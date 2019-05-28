package org.broadinstitute.hellbender.tools.walkers.variantutils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFEncoder;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.IOException;
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
     * alternate_bases, // req ** concat all alt bases with '|' delimiter
     * alternate_bases.alt, // req
     * alternate_bases.AS_RAW_MQ, // req
     * alternate_bases.AS_MQ_DP, // req
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

        position { // Required
             public String getColumnValue(final VariantContext variant) {
                 return String.valueOf(variant.getStart());
            }
        },

        ref { // Required
            public String getColumnValue(final VariantContext variant) {
                final String referenceBase = variant.getReference().getBaseString();
                if (referenceBase == null) {
                    throw new IllegalArgumentException("Cannot be missing required value for reference_bases");
                }
                return referenceBase;
            }
        },

        alt {
            //TODO what if this field is null?
            public String getColumnValue(final VariantContext variant) {
                List<String> outList = new ArrayList<>();
                for(Allele a : variant.getAlternateAlleles()) {
                    outList.add(a.getDisplayString());
                }
                return String.join(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR, outList);
            }
        },

        AS_RAW_MQ { // Required
            public String getColumnValue(final VariantContext variant) {
                String out = getAttribute(variant, GATKVCFConstants.AS_RAW_RMS_MAPPING_QUALITY_KEY, null);
                if (out == null) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_RAW_MQ");
                }
                return out;
            }
        },

        AS_RAW_MQRankSum {
            public String getColumnValue(final VariantContext variant) {
                return getAttribute(variant, GATKVCFConstants.AS_RAW_MAP_QUAL_RANK_SUM_KEY, "");
            }
        },

        AS_QUALapprox { // Required
            public String getColumnValue(final VariantContext variant) {
                //TODO find a constant for "AS_QUALapprox"
                String out = getAttribute(variant, "AS_QUALapprox", null);
                if (out == null) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_QUALapprox");
                }
                return out;
            }
        },

        AS_RAW_ReadPosRankSum {
            public String getColumnValue(final VariantContext variant) {
                return getAttribute(variant, GATKVCFConstants.AS_RAW_READ_POS_RANK_SUM_KEY, "");
            }
        },

        AS_SB_TABLE { // Required
            public String getColumnValue(final VariantContext variant) {
                String out = getAttribute(variant, GATKVCFConstants.AS_SB_TABLE_KEY, null);
                if (out == null) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_SB_TABLE");
                }
                return out;
            }
        },

        AS_VarDP { // Required
            public String getColumnValue(final VariantContext variant) {
                //TODO find a constant for "AS_VarDP"
                String out = getAttribute(variant, "AS_VarDP", null);
                if (out == null) {
                    throw new IllegalArgumentException("Cannot be missing required value for alternate_bases.AS_VarDP");
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
                return variant.getGenotype(0).hasAD() ? Arrays.stream(variant.getGenotype(0).getAD())
                        .mapToObj(String::valueOf)
                        .collect(Collectors.joining(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR)) : "";
            }
        },

        call_DP {
            public String getColumnValue(final VariantContext variant) {
                return variant.getGenotype(0).hasDP() ? String.valueOf(variant.getGenotype(0).getDP()): "";
            }
        },

        call_GQ { // Required
            public String getColumnValue(final VariantContext variant) {
                if (!variant.getGenotype(0).hasGQ()) {
                    throw new IllegalArgumentException("Cannot be missing required value for call.GQ");
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


    public static List<String> createVariantRow(final VariantContext variant) {
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
