package com.idorsia.research.chem.hyperspace3d.model;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.MessageDigest;
import java.util.HexFormat;
import java.util.List;

/** Versioned directed pose-query contract; candidate compatibility is checkpoint-bound. */
public record SharedCandidatePoseBundle(Path directory, JsonNode manifest, String hash) {
    public static SharedCandidatePoseBundle load(Path directory, DeepSpaceModelBundle candidate) {
        try {
            JsonNode m = new ObjectMapper().readTree(directory.resolve("manifest.json").toFile());
            var sums=new java.util.HashMap<String,String>();
            for(String line:Files.readAllLines(directory.resolve("checksums.sha256"))) {
                String[] fields=line.trim().split("\\s+",2);
                require(fields.length==2 && fields[0].matches("[a-f0-9]{64}"),"checksum registry");
                require(sums.put(fields[1],fields[0])==null,"duplicate checksum entry");
            }
            for(String name:List.of("manifest.json","query.onnx","comparator.onnx","feature-schema.json"))
                require(sha256(directory.resolve(name)).equals(sums.get(name)),name+" registry checksum");
            require(m.path("artifactType").asText().equals("deepspace7-shared-candidate-pose-onnx-bundle"), "artifact type");
            require(m.path("architecture").asText().equals("shared_candidate_pose_v1"), "architecture");
            require(m.path("formatVersion").asInt() == 1, "format version");
            require(m.path("graphFeatureContract").asText().equals("deepspace7-v2"), "graph contract");
            require(m.path("distanceUnits").asText().equals("angstrom")
                    && m.path("distanceRepresentation").asText().equals("euclidean_heavy_atom_same_order"), "distance contract");
            require(m.path("directed").asBoolean() && m.path("scoreSemantics").asText().equals("learned_pose_to_ensemble"), "score semantics");
            for (var entry : java.util.Map.of("maxAtoms",32,"atomFeatureDim",50,"pairFeatureDim",36,
                    "embeddingDim",128,"nodeLatentDim",16,"pairLatentDim",14).entrySet())
                require(m.path(entry.getKey()).asInt() == entry.getValue(), entry.getKey());
            require(m.path("availableTargets").equals(new ObjectMapper().valueToTree(
                    List.of("phesa_total","phesa_pharmacophore","phesa_shape"))), "target order");
            require(m.path("supportedPpWeightRange").equals(new ObjectMapper().valueToTree(List.of(0.0,1.0))), "weight range");
            require(m.path("canonicalPredictorSha256").asText().equals(candidate.manifest().predictorCheckpointSha256), "canonical predictor provenance");
            require(m.path("foundationCheckpointSha256").asText().equals(candidate.manifest().foundationCheckpointSha256), "foundation provenance");
            require(m.path("requiredCandidateEncoderSha256").asText().equals(candidate.manifest().encoderSha256), "candidate encoder provenance");
            for (String key : List.of("poseCheckpointSha256","trainingLatentManifestSha256"))
                require(m.path(key).asText().matches("[a-f0-9]{64}"), key);
            for (String name : List.of("query","comparator","feature-schema")) {
                String key = name.equals("feature-schema") ? "featureSchemaSha256" : name+"Sha256";
                require(sha256(directory.resolve(name + (name.equals("feature-schema") ? ".json" : ".onnx")))
                        .equals(m.path(key).asText()), name+" checksum");
            }
            require(m.path("featureSchemaSha256").asText().equals(candidate.manifest().featureSchemaSha256), "feature schema provenance");
            return new SharedCandidatePoseBundle(directory.toAbsolutePath(), m, sha256(directory.resolve("manifest.json")));
        } catch (Exception e) { throw new DeepSpaceInferenceException("Invalid shared-candidate pose bundle: " + e.getMessage(), e); }
    }
    private static void require(boolean condition, String field) {
        if (!condition) throw new IllegalArgumentException("Incompatible pose bundle: " + field);
    }
    public static String sha256(Path path) throws Exception {
        var digest = MessageDigest.getInstance("SHA-256");
        try (var in = Files.newInputStream(path)) {
            byte[] buffer = new byte[1024*1024];
            for (int n; (n=in.read(buffer)) != -1;) digest.update(buffer,0,n);
        }
        return HexFormat.of().formatHex(digest.digest());
    }
    public DeepSpaceModelManifest scoreManifest() {
        var result = new DeepSpaceModelManifest();
        result.availableTargets = List.of("phesa_total","phesa_pharmacophore","phesa_shape");
        return result;
    }
}
