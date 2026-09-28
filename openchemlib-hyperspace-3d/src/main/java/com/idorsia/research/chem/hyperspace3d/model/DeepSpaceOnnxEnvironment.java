package com.idorsia.research.chem.hyperspace3d.model;

import ai.onnxruntime.OrtEnvironment;
import ai.onnxruntime.OrtException;
import ai.onnxruntime.OrtProvider;
import ai.onnxruntime.OrtSession;

public final class DeepSpaceOnnxEnvironment {
    public enum Device { CPU, CUDA }
    private final OrtEnvironment environment;
    private final Device device;
    private final int cudaDeviceId;

    public DeepSpaceOnnxEnvironment(Device device) { this(device, 0); }

    public DeepSpaceOnnxEnvironment(Device device, int cudaDeviceId) {
        this.environment = OrtEnvironment.getEnvironment("hyperspace3d");
        this.device = device;
        this.cudaDeviceId = cudaDeviceId;
        if (device == Device.CUDA && !OrtEnvironment.getAvailableProviders().contains(OrtProvider.CUDA)) {
            throw new DeepSpaceInferenceException(
                    "CUDA was requested but the CUDA execution provider is unavailable; CPU fallback is disabled");
        }
    }

    OrtEnvironment environment() { return environment; }

    OrtSession open(java.nio.file.Path model) {
        try {
            OrtSession.SessionOptions options = new OrtSession.SessionOptions();
            if (device == Device.CUDA) options.addCUDA(cudaDeviceId);
            else options.addCPU(true);
            return environment.createSession(model.toString(), options);
        } catch (OrtException e) {
            throw new DeepSpaceInferenceException("cannot create ONNX session for " + model, e);
        }
    }
}
