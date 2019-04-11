package com.ecn.ferretmvc.model;
import com.ecn.ferretmvc.model.FerretData;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.sql.Connection;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;



public class AnnotThreading implements Runnable{
//	public static URL url;
//	public static URL urlvep;
//	public String [] Varlist;
	//public static volatile String [] test;
	public static volatile Map<Integer, String> StockLineAnnot;
	private int cpt;
	private String Varid;
	private Connection conn;
	 
		public AnnotThreading (int cpt,String Varid,Connection conn) {
			
		 if(AnnotThreading.StockLineAnnot == null)
			 AnnotThreading.StockLineAnnot = new HashMap<>();
		 
		 this.cpt = cpt; 
		 this.Varid = Varid;
		 this.conn = conn;
		}

		@Override
		public void run() {
			 String geneSymbol = null;
			 String geneId = null;
			 String proteinPos = null;
			 String proteinAcc = null;
			 String fxnName = null;
			 String aa1 = null;
			 String aa2 = null;
			 String protein_end = null;
			 String sift_score = null;
			 String sift_score1 = null;
			 String polyphen_score = null;
			 String polyphen_score1 = null;
			 String sift_prediction = null;
			 String sift_prediction1 = null;
			 String polyphen_prediction = null;
			 String polyphen_prediction1 = null;
			 Boolean passee = false;
			 Boolean passee1 = false;
        	 BufferedReader br = null;
        	 BufferedReader brvep = null;
				String regulomeDB_score = null;
				String rgdb_score_signification = null;


         	   System.out.print("Varid threading : " + Varid);

        	
        	
        	try {
            	System.out.println("\tCe qu'il y'a dans varid dans le call run\t" + Varid);
                URL urlLocation = new URL("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_gene.cgi?connect=&rs=" + Varid);
//                String server = "https://rest.ensembl.org";
//    		    String ext = "/vep/human/id/rs" + Varid+ "?content-type=application/json";
//    		    URL urlvep = new URL(server + ext);

    		    br = null;
//    		    brvep = null;
//    		    URLConnection connection = urlvep.openConnection();
//    		    HttpURLConnection httpConnection = (HttpURLConnection)connection;
//       		
//    		    httpConnection.setRequestProperty("Content-Type", "application/json");
    		    	
    		    try{
    		    	
    		    	br = new BufferedReader(new InputStreamReader(urlLocation.openStream()));

                
                String currentString;

                    while ((currentString = br.readLine()) != null  && !currentString.contains("\"mrnaAcc\" : ")) 
                    {
                        if (currentString.contains("\"geneSymbol\"")) {
                        	geneSymbol = currentString.substring(currentString.indexOf(" : \"") + 4, currentString.indexOf("\","));
                        	System.out.print("\tgeneSymbol" + geneSymbol);
                        	 if (geneSymbol.equals("")) {
                             	geneSymbol = ".";
                             }
                    }
                       
                        
                        if (currentString.contains("\"geneId\"")) {
                        	geneId = currentString.substring(currentString.indexOf(" : \"") + 4, currentString.indexOf("\","));
                        	System.out.print("\tgeneId" + geneId);
                        	if (geneId.equals("")) {
                             	geneId = ".";
                             }
                        }
                       
                        
                        
                        if (currentString.contains("\"proteinPos\"")) {
                        	proteinPos = currentString.substring(currentString.indexOf(" : \"") + 4, currentString.indexOf("\","));
                        	System.out.print("\tproteinPos" + proteinPos);
                        	if (proteinPos.equals("")) {
                             	proteinPos = ".";
                             }
                        	
                        	}
                        	
                        
                        if (currentString.contains("\"proteinAcc\"")) {
                        	proteinAcc = currentString.substring(currentString.indexOf(" : \"") + 4, currentString.indexOf("\","));
                        	System.out.print("\tproteinAcc" + proteinAcc);
                        	if (proteinAcc.equals("")) {
                            	proteinAcc = ".";
                            }
                        }
                        
                    
                        if (currentString.contains("\"aaCode\"") && passee1 == false) {
                        	aa1 = currentString.substring(currentString.indexOf(" : \"") + 4, currentString.indexOf("\","));
                        	System.out.print("\t aacode1" + aa1);
                        	if (aa1.equals("")) {
                        		aa1 = ".";
                             }
                        	
                        	passee1 = true;
                        }
                        if (currentString.contains("\"aaCode\"") && passee1 == true) {
                        	
                        	aa2 = currentString.substring(currentString.indexOf(" : \"") + 4, currentString.indexOf("\","));
                        	
                        	System.out.print("\t aacode2" + aa2);
                        	
                        	if (aa2.equals("")) {
                        		aa2 = ".";
                            
                             }
                        
                        }
                        
                       
                        if (currentString.contains("\"fxnName\"") && passee == false) {
                        	
                        	fxnName = currentString.substring(currentString.indexOf(" : \"") + 4, currentString.indexOf("\","));
                        	System.out.print("\tfxnName" + fxnName);
                        	passee=true;
                        	if (fxnName.equals("")) {
                             	fxnName = ".";
                             }
                        	
                        }
                    }
                  //br.close();

// 	           ***********ICI COMMENCE L'EXTRACTION VEP*************
 	           
//                    System.out.println("****** Content of the URL ********");			
//             	   JSONParser parser = new JSONParser();
//             	   String input;
//             	
//
//             		  brvep = new BufferedReader(new InputStreamReader(urlvep.openStream()));
//					while ((input = brvep.readLine()) != null){
//						   //System.out.println("brvepin" + brvep);
//						  //System.out.println("brin" + br);
//						 //System.out.println("input\t" + input);
//						 
//						 JSONArray array = (JSONArray) parser.parse(input);
//						 JSONObject JO = (JSONObject) array.get(0);
//						   JSONArray transcript_consequences = (JSONArray)JO.get("transcript_consequences");
//
//					        for(Object obj: transcript_consequences){
//					     	   JSONObject object = (JSONObject) obj;
//					     	   protein_end = String.valueOf(object.get("protein_end"));
//					     	 if ((!(protein_end.equals("null")) && protein_end.equals(proteinPos)))
//					     	 {
//					     		 System.out.println("Je rentre dans le if vep");
//					     	System.out.println("protein_end: " + protein_end);
//					     	System.out.println(obj);
//					     sift_score1 = String.valueOf(object.get("sift_score"));
//					     
//					     System.out.println("sift_score1: " + sift_score1);
//					     if(!(sift_score1.equals("null"))){
//					    	 sift_score = sift_score1;
//					     System.out.println("sift_score: " + sift_score);
//					     }
//					     if((sift_score1.equals("null"))){
//					        	 sift_score = ".";
//					         System.out.println("sift_scorep: " + sift_score);
//					         }
//					     
//					     polyphen_score1 = String.valueOf(object.get("polyphen_score"));
//					     System.out.println("polyphen_score1: " + polyphen_score1);
//					     
//					     if(!(polyphen_score1.equals("null"))){
//					     polyphen_score = polyphen_score1;
//					     System.out.println("polyphen_score: " + polyphen_score);
//					     
//					     }
//					     if((polyphen_score1.equals("null"))){
//					         polyphen_score = ".";
//					         System.out.println("polyphen_scorep: " + polyphen_score);
//					         }
//					     
//
//					     sift_prediction1 = String.valueOf(object.get("sift_prediction"));
//					     System.out.println("sift_prediction1: " + sift_prediction1);
//					     if (!(sift_prediction1.equals("null")))
//					     {
//					    	 sift_prediction = sift_prediction1;
//					    	 System.out.println("sift_prediction: " + sift_prediction);
//					     }
//					     if ((sift_score1.equals("null")))
//					     {
//					    	 sift_prediction = ".";
//					    	 System.out.println("sift_predictionp: " + sift_prediction);
//					     }
//					     
//					     polyphen_prediction1 = String.valueOf(object.get("polyphen_prediction"));
//					     System.out.println("polyphen_prediction1: " + polyphen_prediction1);
//					     if (!(polyphen_prediction1.equals("null")))
//					     {
//					    	 polyphen_prediction = polyphen_prediction1;
//					    	 System.out.println("polyphen_prediction: " + polyphen_prediction);
//					     }
//					     if ((polyphen_prediction1.equals("null")))
//					     {
//					    	 polyphen_prediction = ".";
//					    	 System.out.println("polyphen_predictionp: " + polyphen_prediction);
//					     }
//					     	 } 
//					         
//					        }
//					        
//					        if(proteinPos.equals("."))
//					        
//					        {
//					    	  System.out.println("On rentre dans un if protend = .");
//					        	sift_score = ".";
//					        	polyphen_score = ".";
//					        	sift_prediction = ".";
//					        	polyphen_prediction = ".";
//					        }
//						 
//					   }
					// Requesting regulomeDB and preparing explanation for the score
					ResultSet res = conn.createStatement()
							.executeQuery("SELECT score FROM scores_table WHERE variant='rs" + Varid + "';");
					if (res.next()) {
						regulomeDB_score = res.getString("score");
						if (regulomeDB_score.equals("1a") || regulomeDB_score.equals("1b")
								|| regulomeDB_score.equals("1c") || regulomeDB_score.equals("1d")
								|| regulomeDB_score.equals("1e") || regulomeDB_score.equals("1f")) {
							rgdb_score_signification = "(likely to affect binding and linked to expression of a gene target)";
						} else if (regulomeDB_score.equals("2a") || regulomeDB_score.equals("2b")
								|| regulomeDB_score.equals("2c")) {
							rgdb_score_signification = "(likely to affect binding)";
						} else if (regulomeDB_score.equals("3a") || regulomeDB_score.equals("3b")) {
							rgdb_score_signification = "(less likely to affect binding)";
						} else if (regulomeDB_score.equals("4") || regulomeDB_score.equals("5")
								|| regulomeDB_score.equals("6")) {
							rgdb_score_signification = "(minimal binding evidence)";
						} else if (regulomeDB_score.equals("7")) {
							rgdb_score_signification = "(no data)";
						}
					} else {
						regulomeDB_score = "not in the base";
						rgdb_score_signification = "(regulomeDB may not be up to date)";
					}

                    AnnotThreading.StockLineAnnot.put(this.cpt,geneSymbol + "\t" + geneId+ "\t" + fxnName + "\t" + proteinPos + "\t" + aa2 + aa1  + "\t" + proteinAcc+ "\t" + sift_score + "\t" + sift_prediction + "\t"
							+ polyphen_score + "\t" + polyphen_prediction+ "\t" + regulomeDB_score + "\t"
							+ rgdb_score_signification);
                    
    			} catch (IOException e) {
    				
					e.printStackTrace();
					System.out.print("e");
					//if (message : tantantan)
				}catch (Exception e2) {
              	     e2.printStackTrace();
              	   System.out.print("e1");
                }
        	} catch (IOException e) {
				e.printStackTrace();
				System.out.print("e2");
            }
		}

		
		
/*public static void Ecriturefile() throws Exception {

	File Freqtest = new File("C:\\Users\\e300264\\Documents\\Projet_de_Stage2018_Ferret\\Fichiers\\FreqTest.txt");
	PrintWriter frqWrite = new PrintWriter(new BufferedWriter(new FileWriter(Freqtest)));
	frqWrite.write("CHROM\tVARIANT\tPOS\tALLELE1\tALLELE2\tNB_CHR\t1KG_A1_FREQ\t1KG_A2_FREQ\tGENENAME\tGENEID\tFUNCTION\tPROTEINPOS\tAACHANGE\tPROTEINACC\tSIFT_SCORE\tSIFT_PREDICTION\tPOLYPHEN_SCORE\tPOLYPHEN_PREDICTION");
	frqWrite.println("");
	frqWrite.println(CallerRunsPolicyDemo.StockLineAnnot[cpt]);
//	frqWrite.write(geneSymbol + "\t" + geneId+ "\t" + fxnName + "\t" + proteinPos + "\t" + aa2 + aa1  + "\t" + proteinAcc + "\t" + sift_score+ "\t" + sift_prediction+ "\t" + polyphen_score+ "\t" + polyphen_prediction);
	frqWrite.println("");
	frqWrite.close();
}*/



	   }

	 
//public static void main(String... args) {
//Runnable runner =new CallerRunsPolicyDemo().new MyRunnable(url,urlvep);
//  
//ThreadPoolExecutor  executor = new ThreadPoolExecutor(2,3,500, TimeUnit.MILLISECONDS,
// 		 new LinkedBlockingQueue<Runnable>(), new ThreadPoolExecutor.CallerRunsPolicy());
//  
//executor.execute(runner);
//	 executor.shutdown();
//		// Wait until all threads are finish
//		while (!executor.isTerminated()) {
//
//		}
//		System.out.println("\nFinished all threads");
//	}


