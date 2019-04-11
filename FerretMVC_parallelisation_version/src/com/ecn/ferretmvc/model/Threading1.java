
package com.ecn.ferretmvc.model;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.*;

import org.broad.tribble.readers.TabixReader;

import com.sun.jimi.core.*;


public class Threading1 implements Runnable {
	
    int j;
    long startTime;
    ArrayList<InputRegion> sortedQueries;
    String webAddress;
    String ftpAddress;
    int queryNumber;
    public static volatile int variantCounter;
    Integer tempInt;
    int[] querySize;
    String fileName;
    //int cpt;
    public static volatile Map<Integer, Map> StockLineVcfall;
   
 
		public Threading1(int j,long startTime,ArrayList<InputRegion> sortedQueries,String webAddress,String ftpAddress,int queryNumber,int variantCounter,Integer tempInt,int[] querySize, String fileName) {
			if(Threading1.StockLineVcfall == null)
				 Threading1.StockLineVcfall = new HashMap<>();
		    this.j = j;
		    //this.cpt = cpt;
		    this.startTime = startTime ;
		    this.sortedQueries = sortedQueries ;
		    this.webAddress = webAddress;
		    this.ftpAddress = ftpAddress;
		    this.queryNumber = queryNumber;
//		    if(Threading1.variantCounter == 0)
//		    	Threading1.variantCounter = variantCounter;
		    this.tempInt = tempInt;
		    this.querySize = querySize;
		    this.fileName = fileName;			
			 			}
		
		
 
		@Override
		public void run() {
			try{
				//File vcfFile = new File(fileName + "_genotypes.vcf");
				
			webAddress = ftpAddress.replace("$", sortedQueries.get(j).getChr());
            System.out.println("QUERY NUMBER ---->" + queryNumber + "\tSORTED QUERIES ---->" + sortedQueries.get(j).getChr() + "\tWEBADRESS ---->" + webAddress);
            //webAddress = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr" + sortedQueries.get(j).getChr() + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz";
            //webAddress = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr" + sortedQueries.get(j).getChr() + ".phase3_shapeit2_mvncall_integrated_v3plus_nounphased.rsID.genotypes.GRCh38_dbSNP_no_SVs.vcf.gz";
            TabixReader tr = new TabixReader(webAddress);
            System.out.println("tr --->" + tr);
           // BufferedWriter vcfBuffWrite = new BufferedWriter(new FileWriter(vcfFile));
            // Get the iterator
            int diff = (int)((sortedQueries.get(j).getEnd() - sortedQueries.get(j).getStart())/2);
            System.out.println("diff---->" + diff);
            int Mid = (sortedQueries.get(j).getStart() + diff);
            System.out.println("Mid---->" + Mid);
            System.out.println("index j ---->" + this.j);
            
           // new tbx(start, mid);
            
           
            
            
            
            startTime = System.nanoTime();
            TabixReader.Iterator iter = tr.query(sortedQueries.get(j).getChr() + ":" + sortedQueries.get(j).getStart() + "-" + sortedQueries.get(j).getEnd());
            System.out.println("iter --->" + iter);
            System.out.println("WHAT'S UP ? --->" + sortedQueries.get(j).getChr() + ":" + sortedQueries.get(j).getStart() + "-" + sortedQueries.get(j).getEnd());
            
            long endTime = System.nanoTime();
            System.out.println("Tabix iterator timeff: " + (endTime - startTime));
            String s;
            Map<Integer, String> StockLineVcfallLocal = new HashMap<>();
            
            int cpt = 0;
            while ((s = iter.next()) != null) {
                Threading1.variantCounter++;
                System.out.println("variantCounter" + Threading1.variantCounter);
                String[] stringSplit = s.split("\t");
                if (stringSplit[6].equals("PASS")) {
                    if (FerretData.isInteger(stringSplit[1])) {
                        tempInt = Integer.parseInt(stringSplit[1]) - sortedQueries.get(j).getStart();
//                        if (tempInt > 0) {
//                            setProgress((int) ((tempInt + querySize[j]) / (double) querySize[queryNumber] * 99));
//                        }
                    }
                    if (stringSplit[2].equals(".")) {
                        stringSplit[2] = "chr" + sortedQueries.get(j).getChr() + "_" + stringSplit[1];
                        String temp = "";
                        for (int i = 0; i < stringSplit.length; i++) {
                        	temp += stringSplit[i] + "\t";
                        	 
                            //vcfBuffWrite.write(stringSplit[i] + "\t");
                            //System.out.println("StringSplit---->" + stringSplit[i] + "\t");
                            //System.out.println("Thr1----->" +StockLineVcfallLocal + "\t");
                        }
                        StockLineVcfallLocal.put(cpt,temp);
                    } else {
                    	StockLineVcfallLocal.put(cpt,s);
                        //vcfBuffWrite.write(s);
                        //System.out.println("s1----->" + s);
                        //System.out.println("Thr1s1----->" +StockLineVcfallLocal + "\t");
                    }
                    //vcfBuffWrite.newLine();
                } else {
                	StockLineVcfallLocal.put(cpt,s);
                	//System.out.println("s2----->" + s);
                    //System.out.println("Thr1s2----->" +StockLineVcfallLocal + "\t");
                    //vcfBuffWrite.write(s);
                    //System.out.println(s);
                    //vcfBuffWrite.newLine();
                }
                cpt++;
                System.out.println("cpt1--->" + cpt);
            }
            Threading1.StockLineVcfall.put(this.j, StockLineVcfallLocal);
            long endEndTime = System.nanoTime();
            System.out.println("Iteration timeff: " + (endEndTime - endTime));
		} catch (IOException e) {
			
			e.printStackTrace();
			System.out.print("e--------------------");
			}
         	
	            }
			
            
            }

