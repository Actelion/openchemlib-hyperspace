package com.idorsia.research.chem.hyperspace.cleaning;

import java.util.Objects;
import java.util.Random;

/**
 * Parameter object that controls how synthon cleaning is orchestrated.
 */
public final class SynthonCleaningOptions {

    private static final String DEFAULT_FRAGMENT_ID_FORMAT = "%s__clean%d";

    private final int assembliesPerSynthon;
    private final boolean reuseFragmentIdWhenSingleResult;
    private final String fragmentIdFormat;
    private final Long randomSeed;

    private SynthonCleaningOptions(Builder builder) {
        this.assembliesPerSynthon = builder.assembliesPerSynthon;
        this.reuseFragmentIdWhenSingleResult = builder.reuseFragmentIdWhenSingleResult;
        this.fragmentIdFormat = builder.fragmentIdFormat;
        this.randomSeed = builder.randomSeed;
    }

    public static Builder builder() {
        return new Builder();
    }

    public int getAssembliesPerSynthon() {
        return assembliesPerSynthon;
    }

    public boolean reuseFragmentIdWhenSingleResult() {
        return reuseFragmentIdWhenSingleResult;
    }

    public String formatFragmentId(String baseId, int duplicateIndex) {
        Objects.requireNonNull(baseId, "baseId");
        return String.format(fragmentIdFormat, baseId, duplicateIndex + 1);
    }

    public Random newRandom() {
        return randomSeed == null ? new Random() : new Random(randomSeed);
    }

    public static final class Builder {
        private int assembliesPerSynthon = 1;
        private boolean reuseFragmentIdWhenSingleResult = true;
        private String fragmentIdFormat = DEFAULT_FRAGMENT_ID_FORMAT;
        private Long randomSeed;

        public Builder assembliesPerSynthon(int assembliesPerSynthon) {
            if (assembliesPerSynthon < 1) {
                throw new IllegalArgumentException("assembliesPerSynthon must be >= 1");
            }
            this.assembliesPerSynthon = assembliesPerSynthon;
            return this;
        }

        public Builder reuseFragmentIdWhenSingleResult(boolean reuse) {
            this.reuseFragmentIdWhenSingleResult = reuse;
            return this;
        }

        /**
         * Sets the format string that is used when a single synthon expands into multiple
         * cleaned synthons. {@link String#format(String, Object...)} will be invoked with
         * the original fragment id and a 1-based duplicate index.
         */
        public Builder fragmentIdFormat(String fragmentIdFormat) {
            Objects.requireNonNull(fragmentIdFormat, "fragmentIdFormat");
            // basic validation
            String test = String.format(fragmentIdFormat, "fragment", 1);
            if (test.isBlank()) {
                throw new IllegalArgumentException("fragmentIdFormat must produce non-empty ids");
            }
            this.fragmentIdFormat = fragmentIdFormat;
            return this;
        }

        public Builder randomSeed(long randomSeed) {
            this.randomSeed = randomSeed;
            return this;
        }

        public SynthonCleaningOptions build() {
            return new SynthonCleaningOptions(this);
        }
    }
}
