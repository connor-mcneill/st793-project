
# Introduction to codes 
	
 We include the necessary codes to obtain the results shown in Section 5 and Appendix D of the Supplementary Material.
	
1. For the application in Section 5.2,  the result can be obtained by running application.R. For Figure 1 in Appendix C.3 of the Supplementary Material, it is derived based on plot\_power.R.
	    	
2. For the  simulation study in Section 5.1 and Appendix D.1--D.2 of the Supplementary Material, the results shown in the tables are derived based on Main\_code.R and Ser.R. To obtain the figures, some revisions to Main\_code.R and Ser.R are needed and the data can be generated from the revised files.  Figure 1 in Section 5.1 and  Figures 2 and 3 in Appendix D of the Supplementary Material are derived based on plot\_main.R. Figure 4 in Appendix D of the Supplementary Material is derived based on plot\_kernel.R. Figure 5 in Appendix D of the Supplementary Material is derived based on the files in the folder D2\_Plot\_power, where the data can be generated from model\_power.R and Ser\_power.R and the figure can be obtained based on the generated data and Plot\_power.R.
	    	
3. For the simulation study in Appendix D.3--D.5 of the Supplementary Material, the results are derived based on the files in the folder Supplementary.
	Specifically, for Appendix D.3  of the Supplementary Material, Table 3  is produced based on  aggregation.R and Ser\_aggregation.R in the folder D\_3;  For Appendix D.4  of the Supplementary Material, the data are generated based on  classical.R and
Ser\_classical.R in the folder D\_4, and Figure 6 is derived based on the generated data and  plot\_classical.R in the folder D\_4;  For Appendix D.5 of the Supplementary Material, Table 5 is derived based on  sub\_p2.R and Ser\_sub\_p2.R in the folder D\_5;
For Appendix D.6 of the Supplementary Material, Table 6 is derived based on   single.R and Ser\_single.R in the folder D\_6.

	
# Workflow 
	   
Since the workflow  of simulation study follows a   similar procedure, the workflow of the main simulation study shown in Section 5.1 and Appendix D.1--D.2 of the Supplementary Material is introduced as an example. The simulation study is based on Main\_code.R and Ser.R, where parallel computing is used. In the simulation study, different settings are considered.  For each setting, the workflow is as follows.

1. Install and load the  R packages $\texttt{glmnet}$ and $\texttt{parallel}$.
2. Specify the setting of the simulation study in Main\_code.R:
   1. Set the model and $\rho$ by the arguments $\texttt{data.gen}$ and $\texttt{rho}$ in the $\texttt{simulation}$ function (lines 371--379).  The default model used is $\texttt{data.bino}$, representing the logistic model. Other models can be applied according to the requirements of the users, such as using $\texttt{data.pois}$ for Poisson model, $\texttt{data.model1}$ for model 1, etc.
	 2. The distribution of $\pmb{z}$ is determined inside the used $\texttt{data.gen}$. The default distribution used in logistic and Poisson models is $\mathcal{N}(\pmb{0},\pmb{I}).$ The default distribution used in other models is the mixture type.
	 3. Set the type of $\pmb{\Sigma}$ by the argument $\verb|no_block|$   in the $\texttt{setting}$ function (line 351).
	 4. Set the value $b$ for $\pmb{\beta}$ by the argument $\texttt{omega2}$ in the $\texttt{setting}$ function (line 351), with $\texttt{omega2}$ representing $b^2$.
	 5. Set the  value of $(n,p)$ in line 350. The default value is $(n,p)=(400,1000)$.
3. Set the number of simulations. It is determined by the length of $\texttt{index}$ (line 363 in Main\_code.R) times the value of $\verb|no_cores|$ (line 11 in Ser.R). The default number of simulations is 1000, with the length of $\texttt{index}$ being 20 and the   value of $\verb|no_cores|$ being 50.	    	
4. Run the file Ser.R.

	
#  Execution time 
	
 For the  application in Section 5.2, it takes around 3.7 hours to complete on the computer with  Intel Core i7-4790  8 cores + 32 GB RAM.   For the simulation study, the codes are run on the server with the processor of 4 x Intel Xeon(R) E7-4890 v2 60 cores +  1536GB RAM. The execution time is provided below for reference.
	
1. Table 1 in Section 5.1.
 Time: It takes around 6.3 hours  for each setting when $(n,p)=(400,1000)$.
2. Figure 5 in Appendix D.2 of the Supplementary  Material.
Time: It takes around 1.4 hour   for each setting.
3. Appendix D.3 of the Supplementary  Material.
Time: It takes around 20 minutes     for each setting.
4. Appendix D.4 of the Supplementary  Material.
Time: It takes around 5 minutes.
5. Appendix D.5 of the Supplementary  Material.
Time: It takes around 5 minutes for each setting.
6. Appendix D.6 of the Supplementary  Material.
Time: It takes around 17  minutes.    	
 
