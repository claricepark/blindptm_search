package edu.scripps.pms.helper;

/**
 * @author  Robin Park
 * @version $Id: PMSConstants.java,v 1.13 2005/10/11 00:12:05 rpark Exp $
 */

public final class PMSConstants
{
    public static final int SEQUEST=1;
    public static final int PEP_PROBE=2;
    public static final int GUTENTAG=3;

    public static final int SCORE_TYPE_SP=1;
    public static final int SCORE_TYPE_XCORR=2;

    //Env. variable of the application
    public static final String APP_HOME="MSP_HOME";

    //Menu Configuration
    public static final String PMS_CONFIG="pmsConfig";
    public static final String PMS_HOME="MSP_HOME";
    public static final String MENU="MENU";
    public static final String ADMIN_MENU="ADMIN_MENU";
    public static final String MAIN_MENU="MAIN_MENU";
    public static final String SUB_MENU="SUB_MENU";
    public static final String NAME="NAME";
    public static final String LINK="LINK";

    //Project
    public static final String SALIVA_PROJECT="Saliva Proteome";

    //DBSource
    public static final int NCBI=1;
    public static final int NCI=2;
    public static final int OTHERS=3; //such as UNIT_PROT, IPI, SGD

    //Sequest related constants
    public static final String SEQUEST_PARAM="sequest.params";

    //DTASelect related constants
    public static final String DTASELECT_PARAM="DTASelect.params";

    public static final String[] MAIN_MENU_LIST=
        {"Home", "Search Programs", "Analysis Programs", "Database", "Tools", "Useful Links", };

    public static final String[] MAIN_MENU_LINK_LIST=
            {"index.jsp", "SearchPrograms.jsp", "AnalysisPrograms.jsp", "Database.jsp", "Tools.jsp", "UsefulLinks.jsp", };

    public static final String SEARCH_DATA_FOLDER = "/search_data/";

}
