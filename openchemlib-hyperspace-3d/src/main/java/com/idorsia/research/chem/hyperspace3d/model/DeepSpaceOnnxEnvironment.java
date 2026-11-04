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
    private final int intraOpThreads;

    public DeepSpaceOnnxEnvironment(Device device) { this(device, 0); }

    public DeepSpaceOnnxEnvironment(Device device, int cudaDeviceId) {
        this(device,cudaDeviceId,0);
    }

    public DeepSpaceOnnxEnvironment(Device device, int cudaDeviceId, int intraOpThreads) {
        if(intraOpThreads<0)throw new IllegalArgumentException("intraOpThreads must be nonnegative (0 selects runtime default)");
        this.intraOpThreads=intraOpThreads;
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
            try(OrtSession.SessionOptions options = new OrtSession.SessionOptions()) {
                if(intraOpThreads>0)options.setIntraOpNumThreads(intraOpThreads);
                if (device == Device.CUDA) {
                    try(var cuda = new ai.onnxruntime.providers.OrtCUDAProviderOptions(cudaDeviceId)) {
                        cuda.add("use_tf32","0");
                        options.addCUDA(cuda);
                    }
                }
                else options.addCPU(true);
                return environment.createSession(model.toString(), options);
            }
        } catch (OrtException e) {
            throw new DeepSpaceInferenceException("cannot create ONNX session for " + model, e);
        }
    }
}
