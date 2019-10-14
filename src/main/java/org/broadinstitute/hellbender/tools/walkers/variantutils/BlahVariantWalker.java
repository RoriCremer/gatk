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

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.nio.file.Path;

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
    private HashMap<String, SimpleXSVWriter> vetWriterCollection = new HashMap<>(23);
    private HashMap<String, SimpleXSVWriter> petWriterCollection = new HashMap<>(23);
    private SimpleXSVWriter sampleMetadataWriter = null;

    private GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;
    private GenomeLocSortedSet coverageLocSortedSet;
    private SimpleInterval previousInterval;
    private String sampleName;
    private String currentContig;
    private List<SimpleInterval> userIntervals;


    @Argument(fullName = "vet-out-path",
            shortName = "VO",
            doc = "Path to the directory where the variants TSVs should be written")
    public GATKPathSpecifier vetOutput = null;

    @Argument(fullName = "pet-out-path",
            shortName = "PO",
            doc = "Path to the directory where the positions expanded TSVs should be written")
    public GATKPathSpecifier petOutput = null;

    @Argument(fullName = "sample-metadata-out-path",
            shortName = "SMO",
            doc = "Path to where the sample metadata TSV should be written")
    public GATKPathSpecifier sampleMetadataOutput = null;

    @Argument(fullName = "ref-block-gq-to-ignore",
            shortName = "IG",
            doc = "Ref Block GQ band to ignore, bands of 10 e.g 0-9 get combined to 0, 20-29 get combined to 20",
            optional = true)
    public String gqStateToIgnore = "SIXTY";

    @Override
    public boolean requiresIntervals() {
        return true;
    }

    @Override
    public void onTraversalStart() {

        final VCFHeader inputVCFHeader = getHeaderForVariants();
        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples());
        if (samples.numberOfSamples() > 1){
            throw new UserException("This tool can only be run on single sample vcfs");
        }
        sampleName = samples.getSample(0);

        final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();

        userIntervals = intervalArgumentCollection.getIntervals(seqDictionary);

        final GenomeLocSortedSet genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));
        coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));

        try {
            List<String> intervalList = userIntervals.stream().map(interval -> interval.toString())
                    .collect(Collectors.toList());
            String intervalListBlob = StringUtils.join(intervalList, ", ");
            String intervalListMd5 = Utils.calcMD5(intervalListBlob);
            List<String> sampleListHeader = BlahSampleListCreation.getHeaders();
            sampleMetadataWriter = new SimpleXSVWriter(sampleMetadataOutput.toPath(), SEPARATOR);
            sampleMetadataWriter.setHeaderLine(sampleListHeader);
            final List<String> TSVLineToCreateSampleMetadata = BlahSampleListCreation.createSampleListRow(
                    sampleName,
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
        if (currentContig != variantChr ) {//if the contig tsvs don't exist yet -- create them
            // TODO should this be pulled out into a helper method?
            try {
                final String petDirectory = petOutput.getURIString();
                final Path petOutputPathByChr = new GATKPathSpecifier(petDirectory + sampleName + "_" + variantChr + FILETYPE).toPath(); // TODO does this need a separator '/'
                List<String> petHeader = BlahPetCreation.getHeaders();
                final SimpleXSVWriter petWriter = new SimpleXSVWriter(petOutputPathByChr, SEPARATOR);
                petWriter.setHeaderLine(petHeader);
                petWriterCollection.put(variantChr, petWriter);
            } catch (final IOException e) {
                throw new UserException("Could not create pet outputs", e);
            }

            try {
                final String vetDirectory = vetOutput.getURIString(); // TODO this is still stripping the '/'
                final Path vetOutputPathByChr = new GATKPathSpecifier(vetDirectory + sampleName + "_" + variantChr + FILETYPE).toPath();
                List<String> vetHeader = BlahVetCreation.getHeaders();
                final SimpleXSVWriter vetWriter = new SimpleXSVWriter(vetOutputPathByChr, SEPARATOR);
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
            final List<String> TSVLineToCreateVet = BlahVetCreation.createVariantRow(variant);

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
                    TSVLinesToCreatePet = BlahPetCreation.createSpanDelRows(start, end, variant, sampleName);
                } else {
                    TSVLinesToCreatePet = BlahPetCreation.createPositionRows(start, end, variant, sampleName);
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
            for (List<String> TSVLineToCreatePet : BlahPetCreation.createMissingTSV(genomeLoc.getStart(), genomeLoc.getEnd(), sampleName)) {
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
