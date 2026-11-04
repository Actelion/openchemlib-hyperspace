package com.idorsia.research.chem.hyperspace.workbench;

/** Stable application errors, independent of HTTP/MCP transports. */
public final class WorkbenchException extends RuntimeException {
    private final String code;
    private final int status;
    public WorkbenchException(String code,int status,String message){super(message);this.code=code;this.status=status;}
    public String code(){return code;}
    public int status(){return status;}
    public static void require(boolean condition,String message){if(!condition)throw new WorkbenchException("invalid_argument",400,message);}
}
