package org.broadinstitute.hellbender.tools.walkers.variantutils;

import com.fasterxml.jackson.annotation.JsonValue;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;

import java.lang.reflect.Executable;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties( // TODO -- how should this be edited?
        summary = "Creates position expanded table",
        oneLineSummary = "Creates position expanded table",
        programGroup = ExampleProgramGroup.class,
        omitFromCommandLine = true
)
public final class BlahPetCreation {

    /**
     * Expected headers for the Variant Table (VET)
     */
    public enum GQStateEnum {
        VARIANT("v"),
        DELETION("*"), // TODO is this s or * ?
        ZERO("0"),
        TEN("1"),
        TWENTY("2"),
        THIRTY("3"),
        FORTY("4"),
        FIFTY("5"),
        SIXTY("6"),
        MISSING("m");

        String value;

        GQStateEnum(String value) {
            this.value = value;
        }

    }

    public static List<String> createTSV(final VariantContext variant) throws Exception {

        // TODO hardcode the dropped band for the GQ

        List<String> row = new ArrayList<>();
        String sampleName = variant.getSampleNamesOrderedByName().get(0);

        if (variant.isVariant()){
             row.add(variant.getStart() + "," + sampleName + "," + GQStateEnum.VARIANT.toString());
            //if variant is variant and has additional positions--must be a deletion: add `*` state
            for (int i = variant.getStart() + 1 ; i < variant.getEnd(); i++){
                row.add(i + "," + sampleName + "," + GQStateEnum.DELETION.value);
            }
        } else {
            int genotypeQual = variant.getGenotype(0).getGQ();  // ok because we only have one sample
            GQStateEnum state = null;

            if (genotypeQual < 10) {
                state = GQStateEnum.ZERO;
            } else if (genotypeQual < 20) {
                state = GQStateEnum.TEN;
            } else if (genotypeQual < 30) {
                state = GQStateEnum.TWENTY;
            } else if (genotypeQual < 40) {
                state = GQStateEnum.THIRTY;
            } else if (genotypeQual < 50) {
                state = GQStateEnum.FORTY;
            } else if (genotypeQual < 60) {
                state = GQStateEnum.FIFTY;
            } else {
                throw new Exception("BLHLBHLBHALBHALBHLBHL");
            }

            for (int i = variant.getStart(); i < variant.getEnd(); i++){ // break up ref blocks
                row.add(i + "," + sampleName + "," + state.value);
            }
        }

        return row;
    }

    public static List<String> createMissingTSV(int start, int end, String sampleName) {
        List<String> row = new ArrayList<>();

        for (int i = start; i < end; i ++){
            row.add(i + "," + sampleName + "," + GQStateEnum.MISSING.value);
        }

        return row;
    }
}
