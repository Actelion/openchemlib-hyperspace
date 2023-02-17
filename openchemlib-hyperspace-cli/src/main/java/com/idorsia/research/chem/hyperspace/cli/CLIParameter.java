package com.idorsia.research.chem.hyperspace.cli;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;

import java.util.function.Consumer;

public class CLIParameter extends Option {

    private final String defaultValue;
    private final Class    type;

    private String parsedString;
    private Object parsedValue;

    /**
     *
     * @param name
     * @param hasParameter
     * @param description
     * @param defaultValue
     */
    public CLIParameter(String name, boolean hasParameter,String description, String defaultValue, Class type) {
        super(name,hasParameter,description);
        this.defaultValue = defaultValue;
        this.type = type;

        if(type != null) {
            if (!(type.equals(String.class) || type.equals(Integer.class) || type.equals(Double.class) || type.equals(Boolean.class) )) {
                throw new Error("Unknown type.. Must be one of String, Integer, Double, Boolean");
            }
        }

    }

    public String getDefaultValue() {
        return this.defaultValue;
    }

    /**
     * Parses config and calls setter function
     * @param cli
     * @return
     */
    public String parseConfig(CommandLine cli) {
        if(this.hasArg()) {
            String val = cli.getOptionValue(this,this.defaultValue);
            this.parsedString = val;
            if(this.type.equals( String.class  )) { parsedValue  = val;}
            if(this.type.equals( Integer.class )) { parsedValue  = Integer.parseInt(val); }
            if(this.type.equals( Double.class  )) { parsedValue  = Double.parseDouble(val);}
            if(this.type.equals( Boolean.class )) { parsedValue  = Boolean.parseBoolean(val);}
            return val;
        }
        else {
            this.parsedValue = true;
        }
        return null;
    }

    public String getStringValue() {return this.parsedString;}
    public int getIntValue() {return (Integer) this.parsedValue;}
    public double getDoubleValue() {return (Double) this.parsedValue;}
    public boolean getBooleanValue() {return (Boolean) this.parsedValue;}

}
