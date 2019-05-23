package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.IOException;
import java.util.*;

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
public final class BlahVariantWalker extends ReferenceWalker {
    static final Logger logger = LogManager.getLogger(BlahVariantWalker.class);

    private final int GQ_CUTOFF = 60;
    private final char SEPARATOR = '\t';
    private SimpleXSVWriter vetWriter = null;
    private SimpleXSVWriter petWriter = null;

    private String sampleName = null;
    private Set<Integer> missingPositions = null;
    private VariantContext currentVariantContext = null;


    @Argument(fullName = "vet-table-out-path",
            shortName = "VO",
            doc="Path to where the variants expanded table should be written")
    public GATKPathSpecifier vetOutput = null;

    @Argument(fullName = "pet-table-out-path",
            shortName = "PO",
            doc="Path to where the positions table should be written")
    public GATKPathSpecifier petOutput = null;

    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            doc="variants to count overlaps of")
    private FeatureInput<VariantContext> variants;


    @Override
    public void onTraversalStart() {

        missingPositions = new HashSet<>();

        try {
            List<String> vetHeader = BlahVetCreation.getHeaders();
            vetWriter = new SimpleXSVWriter(vetOutput.toPath(), SEPARATOR);
            vetWriter.setHeaderLine(vetHeader);

        } catch (final IOException e) {
            throw new IllegalArgumentException("Current variant is missing required fields", e);
        }
        try {
            List<String> petHeader = BlahPetCreation.getHeaders();
            petWriter = new SimpleXSVWriter(petOutput.toPath(), SEPARATOR);
            petWriter.setHeaderLine(petHeader);

        } catch (final IOException e) {
            throw new IllegalArgumentException("Current variant is missing required fields", e);
        }
    }

    @Override
    public void apply(final ReferenceContext referenceContext, final ReadsContext readsContext, final FeatureContext featureContext) {
        if(featureContext.getValues(variants).isEmpty()){
            if (sampleName == null){
                missingPositions.add(referenceContext.getInterval().getStart());
            } else {
                petWriter.getNewLineBuilder().setRow(BlahPetCreation.createMissingTSVRow(referenceContext.getInterval().getStart(), sampleName)).write();
            }
        }
        else {
            VariantContext variant = featureContext.getValues(variants).get(0);

            if (sampleName == null) {
                sampleName = variant.getGenotype(0).getSampleName();
            }

            if (variant.equals(currentVariantContext)) {
                return;
            }
            currentVariantContext = variant;

            // create VET output
            if (!variant.isReferenceBlock()) {
                final List<String> TSVLineToCreateVet = BlahVetCreation.createVariantRow(variant);

                // write the variant to the XSV
                SimpleXSVWriter.LineBuilder vetLine = vetWriter.getNewLineBuilder();
                vetLine.setRow(TSVLineToCreateVet);
                vetLine.write();
            }
            // create PET output
            if (variant.getGenotype(0).getGQ() < 100) {
                List<List<String>> TSVLinesToCreatePet;
                TSVLinesToCreatePet = BlahPetCreation.createPositionRows(variant);

                // write the position to the XSV
                for (List<String> TSVLineToCreatePet : TSVLinesToCreatePet) {
                    petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
                }
            }
        }
    }

    @Override
    public Object onTraversalSuccess(){
        for (int position : missingPositions){
            petWriter.getNewLineBuilder().setRow(BlahPetCreation.createMissingTSVRow(position, sampleName)).write();
        }
        return 0;
    }

    @Override
    public void closeTool() {
        if (vetWriter != null ) {
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
    }
}
