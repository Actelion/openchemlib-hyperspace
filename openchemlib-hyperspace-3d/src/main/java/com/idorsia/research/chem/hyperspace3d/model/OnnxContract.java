package com.idorsia.research.chem.hyperspace3d.model;

import ai.onnxruntime.NodeInfo;
import ai.onnxruntime.OnnxJavaType;
import ai.onnxruntime.OrtException;
import ai.onnxruntime.OrtSession;
import ai.onnxruntime.TensorInfo;
import java.util.Arrays;
import java.util.Map;

final class OnnxContract {
    private OnnxContract() {}

    static void requireInput(OrtSession session, String name, OnnxJavaType type,
                             long... fixedTrailingShape) {
        try {
            require(session.getInputInfo(), name, type, fixedTrailingShape);
        } catch (OrtException e) {
            throw new DeepSpaceInferenceException("cannot inspect ONNX input " + name, e);
        }
    }

    static void requireOutput(OrtSession session, String name, OnnxJavaType type,
                              long... fixedTrailingShape) {
        try {
            require(session.getOutputInfo(), name, type, fixedTrailingShape);
        } catch (OrtException e) {
            throw new DeepSpaceInferenceException("cannot inspect ONNX output " + name, e);
        }
    }

    private static void require(Map<String, NodeInfo> info, String name, OnnxJavaType type,
                                long[] fixedTrailingShape) {
        NodeInfo node = info.get(name);
        if (node == null || !(node.getInfo() instanceof TensorInfo tensor)
                || tensor.type != type) {
            throw new DeepSpaceInferenceException("ONNX tensor " + name + " has the wrong type");
        }
        long[] shape = tensor.getShape();
        if (shape.length != fixedTrailingShape.length + 1 || shape[0] >= 0) {
            throw new DeepSpaceInferenceException(
                    "ONNX tensor " + name + " must have a dynamic batch dimension: "
                            + Arrays.toString(shape));
        }
        for (int i = 0; i < fixedTrailingShape.length; i++) {
            if (shape[i + 1] != fixedTrailingShape[i]) {
                throw new DeepSpaceInferenceException(
                        "ONNX tensor " + name + " has an incompatible shape: "
                                + Arrays.toString(shape));
            }
        }
    }
}
