package com.idorsia.research.chem.hyperspace.io;

import com.idorsia.research.chem.hyperspace.SynthonSpace;
import com.idorsia.research.chem.hyperspace.rawspace.RawSynthonSpace;

/**
 * Imports synthon definitions from vendor TSV files into {@link RawSynthonSpace} dumps.
 */
public class RawSynthonSpaceImporter {

    public Result importEnamine(SynthonSpaceParser3.EnamineOptions options) throws Exception {
        SynthonSpaceParser3.ParsedSpace parsed = SynthonSpaceParser3.parseEnamine(options);
        return new Result(parsed.rawSpace(), parsed.synthonSpace());
    }

    public Result importCsvDirectory(SynthonSpaceParser3.CsvDirectoryOptions options) throws Exception {
        SynthonSpaceParser3.ParsedSpace parsed = SynthonSpaceParser3.parseCsvDirectory(options);
        return new Result(parsed.rawSpace(), parsed.synthonSpace());
    }

    public static final class Result {
        private final RawSynthonSpace rawSpace;
        private final SynthonSpace synthonSpace;

        public Result(RawSynthonSpace rawSpace, SynthonSpace synthonSpace) {
            this.rawSpace = rawSpace;
            this.synthonSpace = synthonSpace;
        }

        public RawSynthonSpace getRawSpace() {
            return rawSpace;
        }

        public SynthonSpace getSynthonSpace() {
            return synthonSpace;
        }
    }
}
