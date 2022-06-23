package edu.scripps.pms.util.dtaselect;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import gnu.trove.TIntDoubleHashMap;

public class ModList {

    private String sequence;
    private TIntDoubleHashMap map = new TIntDoubleHashMap(10);
    private double value;

    public ModList(String sequence, char symbol, double value) {
        this.sequence = sequence;
        this.value = value;

        traverse(sequence, symbol);
    }

    private String traverse(String str, char symbol)
    {
        int index = str.indexOf(symbol);

        if (index > 0) {
            StringBuffer sb = new StringBuffer(str);
            map.put(index+1, value);

            return traverse(sb.deleteCharAt(index).toString(), symbol);
        }

        return str;
    }

    public TIntDoubleHashMap getMap()
    {
        return map;
    }

    public static void main(String args[])
    {
        String str = "SADFSDA#FDSAFDS#FDSA#";

        ModList list = new ModList(str, '#', 80);

        int[] keys=list.getMap().keys();
        System.out.println(keys.length);
        for(int i=0;i<keys.length;i++)
        {
            System.out.println("===><>" + keys[i] + " " + list.getMap().get(keys[i]));
        }

    }
}
