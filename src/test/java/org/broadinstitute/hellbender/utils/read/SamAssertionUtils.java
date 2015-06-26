package org.broadinstitute.hellbender.utils.read;

import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.SamFileValidator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Collection of utilities for making common assertions about SAM files for unit testing purposes.
 */
public final class SamAssertionUtils {

    public static void assertSamsEqual(final File sam1, final File sam2, ValidationStringency validationStringency) throws IOException {
        try (final SamReader reader1 = SamReaderFactory.makeDefault().validationStringency(validationStringency).open(sam1);
             final SamReader reader2 = SamReaderFactory.makeDefault().validationStringency(validationStringency).open(sam2)) {
            final SamComparison comparison = new SamComparison(reader1, reader2);
            final boolean equal = comparison.areEqual();
            Assert.assertTrue(equal, "SAM file output differs from expected output");
        }
    }

    public static void assertSamsEqual(final File sam1, final File sam2) throws IOException {
        assertSamsEqual(sam1, sam2, ValidationStringency.DEFAULT_STRINGENCY);
    }

    public static void assertSamsNonEqual(final File sam1, final File sam2) throws IOException {
        try (final SamReader reader1 = SamReaderFactory.makeDefault().open(sam1);
             final SamReader reader2 = SamReaderFactory.makeDefault().open(sam2)) {
            final SamComparison comparison = new SamComparison(reader1, reader2);
            final boolean equal = comparison.areEqual();
            Assert.assertFalse(equal, "SAM files are expected to differ, but they do not");
        }
    }

    public static void assertSamValid(final File sam) throws IOException {
        try (final SamReader samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(sam)) {
            final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);
            validator.setIgnoreWarnings(true);
            validator.setVerbose(true, 1000);
            validator.setErrorsToIgnore(Arrays.asList(SAMValidationError.Type.MISSING_READ_GROUP));
            final boolean validated = validator.validateSamFileVerbose(samReader, null);
            Assert.assertTrue(validated, "SAM file validation failed");
        }
    }

}
