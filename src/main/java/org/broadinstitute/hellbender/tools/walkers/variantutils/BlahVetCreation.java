package org.broadinstitute.hellbender.tools.walkers.variantutils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class BlahVetCreation {

    // GT:AD:DP:GQ:PL:SB       0/1:232,19,0:251:8:8,0,10988,706,11042,11748:143,89,0,19

    /**
     * Expected headers for the Variant Table (VET)
     * start_position, // req
     * reference_bases, // req
     * alternate_bases, // req ** concat all alt bases with '|' delimiter // TODO--wait, why is this here? looks like its in v1, but not v2
     * alternate_bases.alt, // req
     * alternate_bases.AS_RAW_MQ, // req
     * alternate_bases.AS_MQ_DP, // req // TODO this doesn't exist
     * alternate_bases.AS_RAW_MQRankSum,
     * alternate_bases.AS_QUALapprox, // req
     * alternate_bases.AS_RAW_ReadPosRankSum,
     * alternate_bases.AS_SB_TABLE, // req
     * alternate_bases.AS_VarDP, // req
     * call, // req **
     * call_name, // req
     * call_genotype, // req
     * call_AD,
     * call_DP,
     * call_GQ, // req
     * call_PGT,
     * call_PID,
     * call_PL
     */
    public enum HeaderFieldEnum {
        // TODO is this where the validation step (required vs not) lives  -- fail if there is missing data for a required field
        // and just leave it empty if not required
        // TODO whenever we strip out data, if the data is not as expected, throw an error

        position { // Required
             public String getColumnValue(final VariantContext variant) {
                 return String.valueOf(variant.getStart());
            }
        },

        ref { // Required
            public String getColumnValue(final VariantContext variant) {
                final String referenceBase = variant.getReference().getBaseString();
                if (referenceBase == null) {
                    throw new IllegalArgumentException("Cannot be missing required value for reference_bases"); // TODO, should this be UserException too?
                }
                return referenceBase;
            }
        },

        alt { // remove "<NON_REF>"
            //TODO what if this field is null and if <NON_REF> is not there--throw an error
            public String getColumnValue(final VariantContext variant) {
                List<String> outList = new ArrayList<>();
                for(Allele a : variant.getAlternateAlleles()) {
                    if (!a.isNonRefAllele()) { // TODO unit test this
                        outList.add(a.getDisplayString());
                    }
                }
                return String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, outList);
            }
        },

        AS_RAW_MQ {
            // Required
            // TODO these are string (floats) -- make them ints --- Louis will know where the helper method is?
            // TODO sci notation?
            public String getColumnValue(final VariantContext variant) {
                String out = getAttribute(variant, GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, null);
                if (out == null) {
                    throw new UserException("Cannot be missing required value for alternate_bases.AS_RAW_MQ");
                }
                if (out.endsWith("|0.00")) {
                    out = out.substring(0, out.length() - 5);
                } else {
                    throw new UserException("Expected AS_RAW_MQ value to end in |0.00");
                }
                return out;
            }
        },

        AS_RAW_MQRankSum { // TODO -- maybe rely on 1/1 for call_GT, also get rid of the | at the beginning
            public String getColumnValue(final VariantContext variant) {
                String out =  getAttribute(variant, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY, "");
                if (out.contentEquals("||") || out.contentEquals("|||") ) {
                    out = " "; //  TODO is this better than null?
                    return out;
                }
                if (out.startsWith("|")) {
                    out = out.substring(1);
                } else {
                    throw new UserException("Expected AS_RAW_MQRankSum value to begin with a |");
                }
                if (out.endsWith("|NaN")) {
                    out = out.substring(0, out.length() - 4);
                } else {
                    throw new UserException("Expected AS_RAW_MQRankSum value to be ||, ||| or to end in |NaN");
                }
                return out;

            }
        },

        AS_QUALapprox { // Required
            public String getColumnValue(final VariantContext variant) {
                //TODO find a constant for "AS_QUALapprox"
                String out = getAttribute(variant, "AS_QUALapprox", null);
                if (out == null) {
                    throw new UserException("Cannot be missing required value for alternate_bases.AS_QUALapprox");
                }
                if (out.contentEquals("||") || out.contentEquals("|||") ) {
                    out = " "; //  TODO is this better than null? WAIT WHY IS THIS EVER A ||| VAL IF ITS REQUIRED?
                    return out;
                }
                if (out.startsWith("|")) {
                    out = out.substring(1);
                } else {
                    throw new UserException("Expected AS_RAW_MQRankSum value to begin with a |");
                }
                if (out.endsWith("|0")) {
                    out = out.substring(0, out.length() - 2);
                } else {
                    throw new UserException("Expected AS_QUALapprox value to be ||, ||| or to end in |0");
                }
                return out;
            }
        },

        AS_RAW_ReadPosRankSum {  // TODO -- maybe rely on 1/1 for call_GT
            public String getColumnValue(final VariantContext variant) {
                String out =  getAttribute(variant, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY, "");
                if (out.contentEquals("||") || out.contentEquals("|||") ) {
                    out = " "; // TODO is this better than null?
                    return out;
                }
                if (out.startsWith("|")) {
                    out = out.substring(1);
                } else {
                    throw new UserException("Expected AS_RAW_ReadPosRankSum value to begin with a |");
                }
                if (out.endsWith("|NaN")) {
                    out = out.substring(0, out.length() - 4);
                } else {
                    throw new UserException("Expected AS_RAW_ReadPosRankSum value to be ||, ||| or to end in |NaN");
                }
                return out;
            }
        },

        AS_SB_TABLE { // Required // TODO -- we could remove the 0,0 | if the call_GT is 1/1
            public String getColumnValue(final VariantContext variant) {
                String out = getAttribute(variant, GATKVCFConstants.AS_SB_TABLE_KEY, null);
                if (out == null) {
                    throw new UserException("Cannot be missing required value for alternate_bases.AS_SB_TABLE");
                }
                if (out.endsWith("|0,0")) {
                    out = out.substring(0, out.length() - 4);
                } else {
                    throw new UserException("Expected AS_SB_TABLE value to end in |0,0");
                }
                return out;
            }
        },

        AS_VarDP { // Required
            public String getColumnValue(final VariantContext variant) {
                //TODO find a constant for "AS_VarDP"
                String out = getAttribute(variant, "AS_VarDP", null);
                if (out == null) {
                    throw new UserException("Cannot be missing required value for alternate_bases.AS_VarDP");
                }
                return out;
            }
        },

        sample {
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).getSampleName();
            }
        },

        call_GT {
            public String getColumnValue(final VariantContext variant) {
                ArrayList<Integer> allele_indices = new ArrayList<Integer>();
                for (Allele allele : variant.getGenotype(0).getAlleles()){
                    allele_indices.add(GATKVariantContextUtils.indexOfAllele(variant, allele, true, true, true  ));
                }

                if (allele_indices.size() != 2){
                    throw new IllegalArgumentException("GT doesnt have two alleles");
                }
                return variant.getGenotype(0).isPhased() ? StringUtils.join(allele_indices, VCFConstants.PHASED) : StringUtils.join(allele_indices, VCFConstants.UNPHASED) ;
            }
        },

        call_AD {
            public String getColumnValue(final VariantContext variant) {
                String out = variant.getGenotype(0).hasAD() ? Arrays.stream(variant.getGenotype(0).getAD())
                        .mapToObj(String::valueOf)
                        .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) : "";
                if (out.endsWith(",0")) {
                    out = out.substring(0, out.length() - 2);
                } else {
                    throw new UserException("Expected call_AD to have a final value of 0");
                }
                return out;
            }
        },

        call_DP { // TODO ask Laura if we can drop whole column since it looks similar to AS_VarDP
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).hasDP() ? String.valueOf(variant.getGenotype(0).getDP()): "";
            }
        },

        call_GQ { // Required
            public String getColumnValue(final VariantContext variant) {
                if (!variant.getGenotype(0).hasGQ()) {
                    throw new UserException("Cannot be missing required value for call.GQ");
                }
                return  String.valueOf(variant.getGenotype(0).getGQ());
            }
        },

        CALL_PGT {
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).hasAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY) ? String.valueOf(variant.getGenotype(0).getAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_GT_KEY)) : "";
            }
        },

        call_PID {
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).hasAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY) ? String.valueOf(variant.getGenotype(0).getAnyAttribute(GATKVCFConstants.HAPLOTYPE_CALLER_PHASING_ID_KEY)) : "";
            }
        },

        call_PL {
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).hasPL() ? Arrays.stream(variant.getGenotype(0).getPL())
                        .mapToObj(String::valueOf)
                        .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) : "";
            }
        };

        public String getColumnValue(final VariantContext variant) {
            throw new IllegalArgumentException("Not implemented");
        }

        private static String getAttribute(VariantContext vc, String key, String defaultValue){
            Object attr = vc.getAttribute(key);
            if ( attr == null ) return defaultValue;
            if ( attr instanceof String ) return (String)attr;
            if ( attr instanceof List) return StringUtils.join((List)attr, VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);
            return String.valueOf(attr); // throws an exception if this isn't a string
        }
    }


    public static List<String> createVariantRow(final VariantContext variant) { // TODO  throws UserException ?
        List<String> row = new ArrayList<>();

        for ( final HeaderFieldEnum fieldEnum : HeaderFieldEnum.values() ) {
            row.add(fieldEnum.getColumnValue(variant));
        }
        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(HeaderFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
