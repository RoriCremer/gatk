package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
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
import org.broadinstitute.hellbender.tools.variantdb.BlahVetArrayCreation;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

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
    //private HashMap<String, SimpleXSVWriter> vetWriterCollection = new HashMap<>(26);
    //private HashMap<String, SimpleXSVWriter> petWriterCollection = new HashMap<>(26);
    private SimpleXSVWriter vetWriter = null;
    private SimpleXSVWriter petWriter = null;
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

    // To determine which directory (and ultimately table) the sample's data will go into
    // Since tables have a limited number of samples (default is 4k)
    public int getSampleDirectoryNumber(String sampleId, int sampleMod) { // this is based on sample id
        // sample ids 1-4000 will go in directory 001
        int sampleIdInt = Integer.valueOf(sampleId); // TODO--should sampleId just get refactored as a long?
        // subtract 1 from the sample id to make it 1-index (or do we want to 0-index?) and add 1 to the dir
        int directoryNumber = Math.floorDiv((sampleIdInt - 1), sampleMod) + 1; // TODO omg write some unit tests
        return directoryNumber;
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

    @Argument(fullName = "is-array",
            shortName = "IA",
            doc = "Flag if the input vcf is an array",
            optional = true)
    public Boolean isArray = false;

    @Override
    public boolean requiresIntervals() {
        return true; // TODO -- do I need to check the boolean flag on this?
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
        //  Because BigQuery only supports partitioning based on timestamp or integer,
        // sample names will be remapped into sample_id integers
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

        // Mod the sample directories
        final int sampleMod = 4000; // TODO hardcoded for now--potentially an input param?
        int sampleDirectoryNumber = getSampleDirectoryNumber(sampleId, sampleMod);
        parentDirectory = parentOutputDirectory.toPath(); // TODO do we need this? More efficient way to do this?

        // If this sample set directory doesn't exist yet -- create it
        final String sampleDirectoryName = String.valueOf(sampleDirectoryNumber);
        final Path sampleDirectoryPath = parentDirectory.resolve(sampleDirectoryName);
        final File sampleDirectory = new File(sampleDirectoryPath.toString());
        if (! sampleDirectory.exists()){
            sampleDirectory.mkdir();
        }


        // If the pet directory inside it doesn't exist yet -- create it
        final String petDirectoryName = "pet";
        final Path petDirectoryPath = sampleDirectoryPath.resolve(petDirectoryName);
        final File petDirectory = new File(petDirectoryPath.toString());
        if (! petDirectory.exists()){
            petDirectory.mkdir();
        }
        // If the vet directory inside it doesn't exist yet -- create it
        final String vetDirectoryName = "vet";
        final Path vetDirectoryPath = sampleDirectoryPath.resolve(vetDirectoryName);
        final File vetDirectory = new File(vetDirectoryPath.toString());
        if (! vetDirectory.exists()){
            vetDirectory.mkdir();
        }

        // If the metadata directory inside it doesn't exist yet, create it
        final String metadataDirectoryName = "metadata";
        final Path metadataDirectoryPath = sampleDirectoryPath.resolve(metadataDirectoryName);
        final File sampleMetadataOutputDirectory = new File(metadataDirectoryPath.toString());
        if (! sampleMetadataOutputDirectory.exists()){
            sampleMetadataOutputDirectory.mkdir();
        }

        // if the pet & vet & metadata tsvs don't exist yet -- create them

        try {
            // Create a pet file to go into the pet dir for _this_ sample
            final String petOutputName = sampleName + petDirectoryName + FILETYPE;
            final Path petOutputPath = petDirectoryPath.resolve(petOutputName);
            // write header to it
            List<String> petHeader = BlahPetCreation.getHeaders();
            petWriter = new SimpleXSVWriter(petOutputPath, SEPARATOR);
            petWriter.setHeaderLine(petHeader);
        } catch (final IOException e) {
            throw new UserException("Could not create pet outputs", e);
        }

        try {
            // Create a vet file to go into the pet dir for _this_ sample
            final String vetOutputName = sampleName + vetDirectoryName + FILETYPE;
            final Path vetOutputPath = vetDirectoryPath.resolve(vetOutputName);
            // write header to it
            List<String> vetHeader = isArray ?  BlahVetArrayCreation.getHeaders(): BlahVetCreation.getHeaders();
            vetWriter = new SimpleXSVWriter(vetOutputPath, SEPARATOR);
            vetWriter.setHeaderLine(vetHeader);
        } catch (final IOException e) {
            throw new UserException("Could not create vet outputs", e);
        }

        try {
            // Create a metadata file to go into the metadata dir for _this_ sample
            // TODO--this should just be one file per sample set?
            final String sampleMetadataName = sampleName + metadataDirectoryName+ FILETYPE;
            final Path sampleMetadataOutputPath = metadataDirectoryPath.resolve(sampleMetadataName);
            // write header to it
            List<String> sampleListHeader = BlahSampleListCreation.getHeaders();
            sampleMetadataWriter = new SimpleXSVWriter(sampleMetadataOutputPath, SEPARATOR);
            sampleMetadataWriter.setHeaderLine(sampleListHeader);
            // write values
            List<String> intervalList = userIntervals.stream().map(interval -> interval.toString())
                    .collect(Collectors.toList());
            String intervalListBlob = StringUtils.join(intervalList, ", ");
            String intervalListMd5 = Utils.calcMD5(intervalListBlob);
            final List<String> TSVLineToCreateSampleMetadata = BlahSampleListCreation.createSampleListRow(
                    sampleName,
                    sampleId,
                    intervalListMd5,
                    BlahPetCreation.GQStateEnum.valueOf(gqStateToIgnore));
            sampleMetadataWriter.getNewLineBuilder().setRow(TSVLineToCreateSampleMetadata).write();

        } catch (final IOException e) {
            throw new UserException("Could not create sample metadata outputs", e);
        }


        final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

        userIntervals = intervalArgumentCollection.getIntervals(seqDictionary);

        // To set up the missing positions
        final GenomeLocSortedSet genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));
        coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
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

        // create VET output
        if (!variant.isReferenceBlock()) {
            int start = variant.getStart();
            int end = variant.getEnd();
            // check to see if this is an array
            if(isArray) {
                // check if the array variant is homref 0/0 and if it is then add it to the PET as an unknown state
                ArrayList<Integer> allele_indices = new ArrayList<Integer>();
                for (Allele allele : variant.getGenotype(0).getAlleles()){
                    allele_indices.add(GATKVariantContextUtils.indexOfAllele(variant, allele, true, true, true  ));
                }
                String GT = variant.getGenotype(0).isPhased() ? org.apache.commons.lang.StringUtils.join(allele_indices, VCFConstants.PHASED) : org.apache.commons.lang.StringUtils.join(allele_indices, VCFConstants.UNPHASED) ;
                if (GT.equals("0/0")) { // TODO is this too hard coded? Also note the shortcuts taken for the createArrayPositionRows --- no walking to the end or GQ state other than Unknown
                    List<List<String>> TSVLinesToCreatePet;
                    TSVLinesToCreatePet = BlahPetCreation.createArrayPositionRows(get_location(variantChr, start), get_location(variantChr, end), variant, sampleId);

                    for (List<String> TSVLineToCreatePet : TSVLinesToCreatePet) {
\                        petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
                    }
                } else {
                    final List<String> TSVLineToCreateVet = BlahVetCreation.createVariantRow(
                            get_location(variantChr, start),
                            variant,
                            sampleId
                    );

                    // write the variant to the XSV
                    SimpleXSVWriter.LineBuilder vetLine = vetWriter.getNewLineBuilder();
                    vetLine.setRow(TSVLineToCreateVet);
                    vetLine.write();

                    // also add to PET
                    List<List<String>> TSVLinesToCreatePet;
                    TSVLinesToCreatePet = BlahPetCreation.createPositionRows(get_location(variantChr, start), get_location(variantChr, end), variant, sampleId);

                    // write the position to the XSV
                    for (List<String> TSVLineToCreatePet : TSVLinesToCreatePet) {
                        petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
                    }
                }

                // TODO do I want to return here so there is there's no additional work putting together the genomeloc etc?
                return;
            }
            // else, it must be an exome or genome!
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
        if (isArray) { return 0; }
        for (GenomeLoc genomeLoc : uncoveredIntervals) {
            final String contig = genomeLoc.getContig();
            // write the position to the XSV
            for (List<String> TSVLineToCreatePet : BlahPetCreation.createMissingTSV(
                    get_location(contig, genomeLoc.getStart()),
                    get_location(contig, genomeLoc.getEnd()),
                    sampleId
            )) {
                petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
            }
        }
        return 0;
    }

    @Override
    public void closeTool() {
        if (vetWriter != null) {
            try {
                vetWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VET writer", e);
            }
        }
        if (petWriter != null) {
            try {
                petWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close PET writer", e);
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
