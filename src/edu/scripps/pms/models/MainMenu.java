package edu.scripps.pms.models;

import java.util.*;
import org.apache.ecs.html.A;
import org.apache.ecs.html.LI;
import org.apache.ecs.html.UL;

/**
 * @author  Robin Park
 * @version $Id: MainMenu.java,v 1.4 2005/03/01 19:37:29 rpark Exp $
 */
public class MainMenu extends Menu
{
    private List subMenu = new Vector();

    /**
     * MainMenu
     */
    public MainMenu()
    {
    }

    public boolean add(SubMenu sub)
    {
        return subMenu.add(sub);
    }

    public List getSubMenu()
    {
        return subMenu;
    }

    public void setSubMenu(List subMenu)
    {
        this.subMenu = subMenu;
    }

    public String getSubMenuHtml()
    {
        UL ul = new UL();
        SubMenu sMenu;

        for(Iterator<SubMenu> itr=getSubMenu().iterator(); itr.hasNext(); )
        {
            sMenu = itr.next();

            if( "".equals(sMenu.getHyperlink()) )
                ul.addElement( new LI(sMenu.getLabel()) );
            else
                ul.addElement( new LI(new A(sMenu.getHyperlink(), sMenu.getLabel())) );
        }

        return ul.toString();
    }

}
