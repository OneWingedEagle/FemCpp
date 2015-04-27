package components;


import java.awt.Font;
import java.net.URL;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;


public class Button extends JButton {
	public Font tfFont = new Font("Times New Roman", 1, 13);
	
	public Button(){
	super();	
	setFont(tfFont);
	}
	public Button(String st, int alignment ){
		super();	
		setText(st);
		setHorizontalAlignment(alignment);
		setFont(tfFont);
		}
	public Button(String str){
		super();	
		setText(str);
		setFont(tfFont);
		}
	
}
