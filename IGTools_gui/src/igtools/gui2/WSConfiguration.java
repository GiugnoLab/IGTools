/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.gui2;

import java.io.File;
import java.io.FileWriter;
import java.util.LinkedList;
import java.util.List;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 *
 * @author vbonnici
 */
public class WSConfiguration {
	
	static{
		System.setProperty("javax.xml.parsers.DocumentBuilderFactory", "com.sun.org.apache.xerces.internal.jaxp.DocumentBuilderFactoryImpl");
	}
    
    private Document xmldoc = null;
    private Node wsxmlnode = null;
    
    public WSConfiguration() throws Exception{
        loadWSConfiguration();
        getWSSequences();
    }
    
    private String configFilePath = System.getProperty("user.home") + File.separator + ".igtools" + File.separator + "gui2"  + File.separator + "wsconfig.xml";
    private void loadWSConfiguration() throws Exception{
        File configFile = new File(configFilePath);
        if(!configFile.exists()){
            (new File(configFile.getParent())).mkdirs();
            configFile.createNewFile();
            FileWriter writer = new FileWriter(configFile);
            writer.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n");
            writer.write("<igtools>");
            writer.write("<gui2>");
            writer.write("<wsconfig>");
            writer.write("</wsconfig>");
            writer.write("</gui2>");
            writer.write("</igtools>");
            writer.flush();
            writer.close();
        }
        
        DocumentBuilderFactory docFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = docFactory.newDocumentBuilder();
        xmldoc = docBuilder.parse(new File(configFilePath));
        
        wsxmlnode = xmldoc.getElementsByTagName("igtools").item(0);
        wsxmlnode = getFirstElement(wsxmlnode, "gui2");
        wsxmlnode = getFirstElement(wsxmlnode, "wsconfig");
        
    }
    
    
    private Node getFirstElement(Node node, String tagname){
        NodeList nodeList = node.getChildNodes();
        for(int i=0; i<nodeList.getLength(); i++){
          Node childNode = nodeList.item(i);
          if(childNode.getNodeType() == Node.ELEMENT_NODE &&
                  childNode.getNodeName().compareTo(tagname) == 0
                  ){
              return childNode;
          }
        }
        return null;
    }
    
    
    
    
    
    public List<WSSequence> getWSSequences(){
        List<WSSequence> seqs = new LinkedList<WSSequence>();
        
        Node seqsContainer = getFirstElement(wsxmlnode, "wssequences");
        if(seqsContainer != null){
            NodeList nodeList = seqsContainer.getChildNodes();
            for(int i=0; i<nodeList.getLength(); i++){
              Node childNode = nodeList.item(i);
              if(childNode.getNodeType() == Node.ELEMENT_NODE &&
                      childNode.getNodeName().compareTo("sequence") == 0
                      ){
                 seqs.add(getWSSequence(childNode));
              }
            }
        }
        
        return seqs;
    }
    private WSSequence getWSSequence(Node node){
        WSSequence seq = new WSSequence();
        
        NodeList nodeList = node.getChildNodes();
        for(int i=0; i<nodeList.getLength(); i++){
          Node childNode = nodeList.item(i);
          if(childNode.getNodeType() == Node.ELEMENT_NODE){
            if(childNode.getNodeName().compareTo("name") == 0 ){
               seq.setName(childNode.getTextContent().trim());
            }
            else if(childNode.getNodeName().compareTo("fastaFile") == 0 ){
               seq.setFASTAFile(childNode.getTextContent().trim());
            }
            else if(childNode.getNodeName().compareTo("b3seqFile") == 0 ){
               seq.setB3seqFile(childNode.getTextContent().trim());
            }
            else if(childNode.getNodeName().compareTo("nelsaFile") == 0 ){
               seq.setNELSAFile(childNode.getTextContent().trim());
            }
            else{}
          }
        }
        return seq;
    }
    
    
    
    
    public void save(List<WSSequence> seqs) throws Exception{
        File configFile = new File(configFilePath);
        if(!configFile.exists()){
            (new File(configFile.getParent())).mkdirs();
            configFile.createNewFile();
        }
       
        FileWriter writer = new FileWriter(configFile);
        writer.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n");
        writer.write("<igtools>\n");
        writer.write("<gui2>\n");
        writer.write("<wsconfig>\n");
        writer.write("<wssequences>\n");


        System.out.println("writing sequences...");
        for(WSSequence seq : seqs){
            writer.write("<sequence>\n");
            writer.write("<name>"+escapeXml(seq.getName())+"</name>\n");
            if(seq.getFASTAFile()!=null)
            writer.write("<fastaFile>"+escapeXml(seq.getFASTAFile())+"</fastaFile>\n");
            if(seq.getB3seqFile()!=null)
            writer.write("<b3seqFile>"+escapeXml(seq.getB3seqFile())+"</b3seqFile>\n");
            if(seq.getNELSAFile()!=null)
            writer.write("<nelsaFile>"+escapeXml(seq.getNELSAFile())+"</nelsaFile>\n");
            writer.write("</sequence>\n");
        }

        writer.write("</wssequences>\n");
        writer.write("</wsconfig>\n");
        writer.write("</gui2>\n");
        writer.write("</igtools>\n");
        writer.flush();
        writer.close();
    }
    
    public String escapeXml(String s) {
        return s.replaceAll("&", "&amp;").replaceAll(">", "&gt;").replaceAll("<", "&lt;").replaceAll("\"", "&quot;").replaceAll("'", "&apos;");
    }
}
