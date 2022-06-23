package edu.scripps.pms.models;


/**
 * @author  Robin Park
 * @version $Id: Menu.java,v 1.3 2005/02/15 18:00:44 rpark Exp $
 */
public class Menu
{
    private boolean visible;
    private String label;
    private String hyperlink;
    private String target;

    public boolean isVisible()
    {
        return visible;
    }

    public String getLabel()
    {
        return label;
    }

    public String getHyperlink()
    {
        return hyperlink;
    }

    public String getTarget()
    {
        return target;
    }

    public void setVisible(boolean visible)
    {
        this.visible = visible;
    }

    public void setLabel(String label)
    {
        this.label = label;
    }

    public void setHyperlink(String hyperlink)
    {
        this.hyperlink = hyperlink;
    }

    public void setTarget(String target)
    {
        this.target = target;
    }

    public String toString()
    {
        return getLabel();
    }

}
