package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.nio.file.Path;
import java.io.File;

/**
 * Example/toy program that shows how to implement the VariantWalker interface. Prints supplied variants
 * along with overlapping reads/reference bases/variants (if present).
 */
@CommandLineProgramProperties(
        summary = "Example tool that prints variants supplied to the specified output file (stdout if none provided), along with overlapping reads/reference bases/variants (if provided)",
        oneLineSummary = "Example tool that prints variants with optional contextual data",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class BlahVariantWalker extends VariantWalker {
    static final Logger logger = LogManager.getLogger(BlahVariantWalker.class);

    private final char SEPARATOR = '\t';
    private final String FILETYPE = ".tsv";
    private HashMap<String, SimpleXSVWriter> vetWriterCollection = new HashMap<>(26);
    private HashMap<String, SimpleXSVWriter> petWriterCollection = new HashMap<>(26);
    private SimpleXSVWriter sampleMetadataWriter = null;

    private GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;
    private GenomeLocSortedSet coverageLocSortedSet;
    private SimpleInterval previousInterval;
    private String sampleName;
    private String sampleId;
    private String currentContig;
    private List<SimpleInterval> userIntervals;

    public long chromAdjustment = 1000000000000L;

    public enum ChromosomeEnum {
        CHR1(1),
        CHR2(2),
        CHR3(3),
        CHR4(4),
        CHR5(5),
        CHR6(6),
        CHR7(7),
        CHR8(8),
        CHR9(9),
        CHR10(10),
        CHR11(11),
        CHR12(12),
        CHR13(13),
        CHR14(14),
        CHR15(15),
        CHR16(16),
        CHR17(17),
        CHR18(18),
        CHR19(19),
        CHR20(20),
        CHR21(21),
        CHR22(22),
        CHRX(23),
        CHRY(24),
        CHRM(25);

        int index;

        ChromosomeEnum(int index) {
            this.index = index;
        }
    }

    public long get_location(String chrom, int position) {
        int chromosomeIndex = ChromosomeEnum.valueOf(chrom.toUpperCase()).index;
        long adjustedLocation = Long.valueOf(chromosomeIndex) * chromAdjustment + Long.valueOf(position);
        return adjustedLocation;
    }

    // Inside the parent directory, a directory for each chromosome will be created, with a pet directory and vet directory in each one.
    // Each pet and vet directory will hold all of the pet and vet tsvs for each sample
    // A metadata directory will be created, with a metadata tsv for each sample

    @Argument(fullName = "vet-pet-out-path",
            shortName = "VPO",
            doc = "Path to the directory where the variants TSVs and positions expanded TSVs should be written")
    public GATKPathSpecifier parentOutputDirectory = null;
    public Path parentDirectory = null;

    @Argument(fullName = "ref-block-gq-to-ignore",
            shortName = "IG",
            doc = "Ref Block GQ band to ignore, bands of 10 e.g 0-9 get combined to 0, 20-29 get combined to 20",
            optional = true)
    public String gqStateToIgnore = "SIXTY";

    @Argument(fullName = "sample-name-mapping",
            shortName = "SNM",
            doc = "Sample name to sample id mapping")
    public File sampleMap;

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {

        // Get sample name
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples());
        if (samples.numberOfSamples() > 1){
            throw new UserException("This tool can only be run on single sample vcfs");
        }
        sampleName = samples.getSample(0);
        try {
            BufferedReader br = new BufferedReader(new FileReader(sampleMap));

            String line; // Reading header, Ignoring
            while ((line = br.readLine()) != null && !line.isEmpty()) {
                String[] fields = line.split(",");
                String name = fields[1];
                if (sampleName.equals(name)) {
                    sampleId = fields[0];
                    break;
                }
            }
            br.close();
            if (sampleId == null) {
                // sampleName not found
                throw new UserException("Sample " + sampleName + " could not be found in sample mapping file");
            }
        } catch (final IOException ioe) { // FileNotFoundException e,
            throw new UserException("Could not find sample mapping file");
        }

        // If the metadata directory doesn't exist yet, create it
        parentDirectory = parentOutputDirectory.toPath();
        final String metadataDirectoryName = "metadata";
        final Path metadataDirectoryPath = parentDirectory.resolve(metadataDirectoryName);
        final File sampleMetadataOutputDirectory = new File(metadataDirectoryPath.toString());

        if (! sampleMetadataOutputDirectory.exists()){
            sampleMetadataOutputDirectory.mkdir();
        }

        // Create a metadata file to go into it for _this_ sample
        final String sampleMetadataName = sampleName + FILETYPE;
        final Path sampleMetadataOutput = metadataDirectoryPath.resolve(sampleMetadataName);


        final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

        userIntervals = intervalArgumentCollection.getIntervals(seqDictionary);

        // To set up the missing positions
        final GenomeLocSortedSet genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));
        coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));

        try {
            List<String> intervalList = userIntervals.stream().map(interval -> interval.toString())
                    .collect(Collectors.toList());
            String intervalListBlob = StringUtils.join(intervalList, ", ");
            String intervalListMd5 = Utils.calcMD5(intervalListBlob);
            List<String> sampleListHeader = BlahSampleListCreation.getHeaders();
            sampleMetadataWriter = new SimpleXSVWriter(sampleMetadataOutput, SEPARATOR);
            sampleMetadataWriter.setHeaderLine(sampleListHeader);
            final List<String> TSVLineToCreateSampleMetadata = BlahSampleListCreation.createSampleListRow(
                    sampleName,
                    sampleId,
                    intervalListMd5,
                    BlahPetCreation.GQStateEnum.valueOf(gqStateToIgnore));
            sampleMetadataWriter.getNewLineBuilder().setRow(TSVLineToCreateSampleMetadata).write();

        } catch (final IOException e) {
            throw new UserException("Could not create sample metadata output", e);
        }
    }

    public void setCoveredInterval( String variantChr, int start, int end) {
        // add interval to "covered" intervals
        // GenomeLocSortedSet will automatically merge intervals that are overlapping when setting `mergeIfIntervalOverlaps`
        // to true.  In a GVCF most blocks are adjacent to each other so they wouldn't normally get merged.  We check
        // if the current record is adjacent to the previous record and "overlap" them if they are so our set is as
        // small as possible while still containing the same bases.
        final SimpleInterval variantInterval = new SimpleInterval(variantChr, start, end);
        final int intervalStart = (previousInterval != null && previousInterval.overlapsWithMargin(variantInterval, 1)) ?
                previousInterval.getStart() : variantInterval.getStart();
        final int intervalEnd = (previousInterval != null && previousInterval.overlapsWithMargin(variantInterval, 1)) ?
                Math.max(previousInterval.getEnd(), variantInterval.getEnd()) : variantInterval.getEnd();

        final GenomeLoc possiblyMergedGenomeLoc = coverageLocSortedSet.getGenomeLocParser().createGenomeLoc(variantInterval.getContig(), intervalStart, intervalEnd);
        coverageLocSortedSet.add(possiblyMergedGenomeLoc, true);
        previousInterval = new SimpleInterval(possiblyMergedGenomeLoc);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        // get the intervals this variant covers
        final GenomeLoc variantGenomeLoc = intervalArgumentGenomeLocSortedSet.getGenomeLocParser().createGenomeLoc(variant.getContig(), variant.getStart(), variant.getEnd());
        final List<GenomeLoc> intervalsToWrite = intervalArgumentGenomeLocSortedSet.getOverlapping(variantGenomeLoc);

        if (intervalsToWrite.size() == 0){
            throw new IllegalStateException("There are no intervals being covered by this variant, something went wrong with interval parsing");
        }

        // take the first interval(assuming this is returned in order) and make sure if its a variant, that it starts at/after the interval start
        // we are going to ignore any deletions that start before an interval.
        if (!variant.isReferenceBlock() && intervalsToWrite.get(0).getStart() > variant.getStart()){
            return;
        }

        // if the only alt allele for a variant is `*`, we ignore it
        if (!variant.isReferenceBlock() &&  variant.getAlternateAlleles().size() == 2 && variant.hasAlternateAllele(Allele.SPAN_DEL)){
            return;
        }

        final String variantChr = variant.getContig();

        // If this contig directory don't exist yet -- create it
        final String contigDirectoryName = variantChr;
        final Path contigDirectoryPath = parentDirectory.resolve(contigDirectoryName);
        final File contigDirectory = new File(contigDirectoryPath.toString());
        if (! contigDirectory.exists()){
            contigDirectory.mkdir();
        }
        // If the pet directory inside it doesn't exist yet -- create it
        final String petDirectoryName = "pet";
        final Path petDirectoryPath = contigDirectoryPath.resolve(petDirectoryName);
        final File petDirectory = new File(petDirectoryPath.toString());
        if (! petDirectory.exists()){
            petDirectory.mkdir();
        }
        // If the vet directory inside it doesn't exist yet -- create it
        final String vetDirectoryName = "vet";
        final Path vetDirectoryPath = contigDirectoryPath.resolve(vetDirectoryName);
        final File vetDirectory = new File(vetDirectoryPath.toString());
        if (! vetDirectory.exists()){
            vetDirectory.mkdir();
        }


        if (currentContig != variantChr ) {// if the pet & vet tsvs don't exist yet -- create them
            try {
                // Create a pet file to go into the pet dir for _this_ sample
                final String petOutputName = sampleName + FILETYPE;
                final Path petOutputPath = petDirectoryPath.resolve(petOutputName);
                // Write to it
                List<String> petHeader = BlahPetCreation.getHeaders();
                final SimpleXSVWriter petWriter = new SimpleXSVWriter(petOutputPath, SEPARATOR);
                petWriter.setHeaderLine(petHeader);
                petWriterCollection.put(variantChr, petWriter);
            } catch (final IOException e) {
                throw new UserException("Could not create pet outputs", e);
            }

            try {
                // Create a vet file to go into the pet dir for _this_ sample
                final String vetOutputName = sampleName + FILETYPE;
                final Path vetOutputPath = vetDirectoryPath.resolve(vetOutputName);
                // Write to it
                List<String> vetHeader = BlahVetCreation.getHeaders();
                final SimpleXSVWriter vetWriter = new SimpleXSVWriter(vetOutputPath, SEPARATOR);
                vetWriter.setHeaderLine(vetHeader);
                vetWriterCollection.put(variantChr, vetWriter);
            } catch (final IOException e) {
                throw new UserException("Could not create vet outputs", e);
            }
            currentContig = variantChr;
        }

        final SimpleXSVWriter vetWriter = vetWriterCollection.get(variantChr);

        // create VET output
        if (!variant.isReferenceBlock()) {
            int start = variant.getStart();
            final List<String> TSVLineToCreateVet = BlahVetCreation.createVariantRow(
                    get_location(variantChr, start),
                    variant,
                    sampleId
            );

            // write the variant to the XSV
            SimpleXSVWriter.LineBuilder vetLine = vetWriter.getNewLineBuilder();
            vetLine.setRow(TSVLineToCreateVet);
            vetLine.write();
        }
        boolean firstInterval = true;
        for (GenomeLoc genomeLoc : intervalsToWrite) {

            int start = Math.max(genomeLoc.getStart(), variant.getStart());
            int end = Math.min(genomeLoc.getEnd(), variant.getEnd());

            // TODO throw an error if start and end are the same?

            // for each of the reference blocks with the GQ to discard, keep track of the positions for the missing insertions
            if (BlahPetCreation.getGQStateEnum(variant.getGenotype(0).getGQ()).equals(BlahPetCreation.GQStateEnum.valueOf(gqStateToIgnore))) {
                // add interval to "covered" intervals
                setCoveredInterval(variantChr, start, end);
            }

            // create PET output if the reference block's GQ is not the one to discard or its a variant
            if (!variant.isReferenceBlock() || !BlahPetCreation.getGQStateEnum(variant.getGenotype(0).getGQ()).equals(BlahPetCreation.GQStateEnum.valueOf(gqStateToIgnore))) {

                // add interval to "covered" intervals
                setCoveredInterval(variantChr, start, end);

                List<List<String>> TSVLinesToCreatePet;
                // handle deletions that span across multiple intervals
                if (!firstInterval && !variant.isReferenceBlock()) {
                    TSVLinesToCreatePet = BlahPetCreation.createSpanDelRows(
                            get_location(variantChr, start),
                            get_location(variantChr, end),
                            variant,
                            sampleId
                    );
                } else {
                    TSVLinesToCreatePet = BlahPetCreation.createPositionRows(
                            get_location(variantChr, start),
                            get_location(variantChr, end),
                            variant,
                            sampleId
                    );
                }

                // write the position to the XSV
                for (List<String> TSVLineToCreatePet : TSVLinesToCreatePet) {
                    final SimpleXSVWriter petWriter = petWriterCollection.get(variantChr);
                    petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
                }
            }
            firstInterval = false;
        }
    }

    @Override
    public Object onTraversalSuccess() {
        final GenomeLocSortedSet uncoveredIntervals = intervalArgumentGenomeLocSortedSet.subtractRegions(coverageLocSortedSet);
        logger.info("MISSING_GREP_HERE:" + uncoveredIntervals.coveredSize());
        logger.info("MISSING_PERCENTAGE_GREP_HERE:" + (1.0 * uncoveredIntervals.coveredSize()) / intervalArgumentGenomeLocSortedSet.coveredSize());
        for (GenomeLoc genomeLoc : uncoveredIntervals) {
            final String contig = genomeLoc.getContig();
            // write the position to the XSV
            for (List<String> TSVLineToCreatePet : BlahPetCreation.createMissingTSV(
                    get_location(contig, genomeLoc.getStart()),
                    get_location(contig, genomeLoc.getEnd()),
                    sampleId
            )) {
                petWriterCollection.get(contig).getNewLineBuilder().setRow(TSVLineToCreatePet).write();
            }
        }
        return 0;
    }

    @Override
    public void closeTool() {
        for (SimpleXSVWriter writer:vetWriterCollection.values()) {
            if (writer != null) {
                try {
                    writer.close();
                } catch (final Exception e) {
                    throw new IllegalArgumentException("Couldn't close VET writer", e);
                }
            }
        }
        for (SimpleXSVWriter writer:petWriterCollection.values()) {
            if (writer != null) {
                try {
                    writer.close();
                } catch (final Exception e) {
                    throw new IllegalArgumentException("Couldn't close PET writer", e);
                }
            }
        }
        if (sampleMetadataWriter != null) {
            try {
                sampleMetadataWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close Sample Metadata writer", e);
            }
        }
    }
}
