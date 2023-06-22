# Data Sets from lýteratures
#
# 1.
#   Zelterman, Daniel. (2004). 
#   Discrete Distributions - Applications in the Health Sciences.
#   John Wýley & Sons, West Sussex.
# 2.
#   Chung, C. S. and Brown, K. S. (1970).
#   Family Studies of Early Childhood Deafness Ascertained through the Clarke School for the Deaf.
#   The American Society of Human Genetics, 22, 630-644.
# 3.
#   Pala F. S, Alkaya F., Tabakcioglu K., Tokatli F., Uzal C., Parlar Þ., Algunes C. (2008).
#   The Effects of Micronuclei with Whole Chromosome on Biological Dose Estimation.
#   Turkish Journal of Biology, 32, 283-290.
# 4.
#   Elwood M. and Coldman A. (1981). 
#   Age of mothers with breast cancer and sex of their children.
#   British Medical Journal, 282, 734.
# 5. 
#   Rao C. R. (1997).
#   Statistics and truth: putting chance to work. 
#   World Scientific Publishing Co., Singapore.
#
# 6. 
#  Hubert Lacey J., Gowers S. G. and Bhat A. V. (1991).
#  Bulimia Nervosa: Family Size, Sibling Sex and Birth Order 
#  A Catchment-Area Study.
#  British Journal of Psychiatry, 158, 491-494.
#
# 7. Haukka, J. K., Suvisaari J. and Lönnqvist J. (2004).
# Family structure and risk factors for schizophrenia: case-sibling study.
# BMC Psychiatry , 4:41, doi:10.1186/1471-244X-4-41.
#
# 8. Omariba D. W. R., Rajulton F., Beaujot R. (2008).
# Correlated mortality of siblings in Kenya: The role of state dependence.
# Demographic Research, 18, 311-336.
#
# 9.Faddy, M. J., and Bosch, R. J. (1999).
# Likelihood-based modeling and analysis of data underdispersed relative to the Poisson distribution. 
# Biometrics, 57, 620-624.



 


#
# Zelterman (2004),pg.144, Tabel 6.1
#
# n	f_n	m_n	0	1	2	3	4	5	6
#
# 1	48	12	36	12	0	0	0	0	0
# 2	23	9	15	7	1	0	0	0	0
# 3	17	19	5	7	3	2	0	0	0
# 4	7	5	3	3	1	0	0	0	0
# 6	5	15	1	0	1	1	1	0	1

tab6.1=matrix(0, nrow=5, ncol=10)
colnames(tab6.1)=c("n", "f_n", "m_n", as.character(0:6))
tab6.1[, 1]=c(1, 2, 3, 4, 6)
tab6.1[, 2]=c(48, 23, 17, 7, 5)
tab6.1[, 3]=c(12, 9, 19, 5, 15)
tab6.1[, 4]=c(36, 15, 5, 3, 1)
tab6.1[, 5]=c(12, 7, 7, 3, 0)
tab6.1[, 6]=c(0, 1, 3, 1, 1)
tab6.1[, 7]=c(0, 0, 2, 0, 1)
tab6.1[, 8]=c(0, 0, 0, 0, 1)
tab6.1[, 9]=c(0, 0, 0, 0, 0)
tab6.1[, 10]=c(0, 0, 0, 0, 1)	


# Zelterman (2004),pg. 144, Tabel 6.2
#
# n	f_n	m_n	0	1	2	3	4+
#
# 0	1	0	1	0	0	0	0
# 1	2	0	2	0	0	0	0
# 2	8	10	1	4	3	0	0
# 3	4	7	0	2	1	1	0
# 4	2	2	0	2	0	0	0
# 5	1	1	0	1	0	0	0
# 6	1	1	0	1	0	0	0
# 7	1	1	0	1	0	0	0
# 8	3	6	0	1	1	1	0
# 9	1	3	0	0	0	1	0


tab6.2=matrix(0, nrow=10, ncol=13)
colnames(tab6.2)=c("n", "f_n", "m_n", as.character(0:9))
tab6.2[, 1]=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
tab6.2[, 2]=c(1, 2, 8, 4, 2, 1, 1, 1, 3, 1)
tab6.2[, 3]=c(0, 0, 10, 7, 2, 1, 1, 1, 6, 3)
tab6.2[, 4]=c(1, 2, 1, 0, 0, 0, 0, 0, 0, 0)
tab6.2[, 5]=c(0, 0, 4, 2, 2, 1, 1, 1, 1, 0)
tab6.2[, 6]=c(0, 0, 3, 1, 0, 0, 0, 0, 1, 0)
tab6.2[, 7]=c(0, 0, 0, 1, 0, 0, 0, 0, 1, 1)
tab6.2[, 8]=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


# Zelterman (2004), pg. 145, Tabel 6.3
#
# n	f_n	m_n	0	1	2	3	4	5	6	7+ 	
#
# 1	267	12	255	12	0	0	0	0	0	0
# 2	285	48	239	44	2	0	0	0	0	0
# 3	202	80	143	41	15	3	0	0	0	0
# 4	110	54	69	30	9	2	0	0	0	0
# 5	104	103	43	34	15	9	3	0	0	0
# 6	50	67	15	18	8	5	3	0	1	0
# 7	21	38	4	4	7	4	2	0	0	0
# 8	12	28	1	2	4	3	1	1	0	0
#


tab6.3=matrix(0, nrow=8, ncol=12)
colnames(tab6.3)=c("n", "f_n", "m_n", as.character(0:8))
tab6.3[, 1]=c(1, 2, 3, 4, 5, 6, 7, 8)
tab6.3[, 2]=c(267, 285, 202, 110, 104, 50, 21, 12)
tab6.3[, 3]=c(12, 48, 80, 54, 103, 67, 38, 28)
tab6.3[, 4]=c(255, 239, 143, 69, 43, 15, 4, 1)
tab6.3[, 5]=c(12, 44, 41, 30, 34, 18, 4, 2)
tab6.3[, 6]=c(0, 2, 15, 9 , 15, 8, 7, 4)
tab6.3[, 7]=c(0, 0, 3, 2, 9, 5, 4, 3)
tab6.3[, 8]=c(0, 0, 0, 0, 3, 3, 2, 1)
tab6.3[, 9]=c(0, 0, 0, 0, 0, 0, 0, 1)
tab6.3[, 10]=c(0, 0, 0, 0, 0, 1, 0, 0)
tab6.3[, 11]=c(0, 0, 0, 0, 0, 0, 0, 0)




# Zelterman (2004), pg.146, Tabel 6.4
#
# n	f_n	m_n	0	1	2	3	4	5	6	7+
#
# 1	1	0	1	0	0	0	0	0	0	0
# 2	8	8	2	4	2	0	0	0	0	0
# 3	13	11	5	5	3	0	0	0	0	0
# 4	5	5	2	1	2	0	0	0	0	0
# 5	5	9	2	1	0	1	0	1	0	0
# 6	4	17	0	0	1	0	1	1	1	0
# 7	1	1	0	1	0	0	0	0	0	0
# 8	1	3	0	0	0	1	0	0	0	0
# 9	1	1	0	1	0	0	0	0	0	0	
# 13	1	2	0	0	1	0	0	0	0	0


tab6.4=matrix(0, nrow=10, ncol=17)
colnames(tab6.4)=c("n", "f_n", "m_n", as.character(0:13))
tab6.4[, 1]=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 13)
tab6.4[, 2]=c(1, 8, 13, 5, 5, 4, 1, 1, 1, 1)
tab6.4[, 3]=c(0, 8, 11, 5, 9, 17, 1, 3, 1, 2)
tab6.4[, 4]=c(1, 2, 5, 2, 2, 0, 0, 0, 0, 0)
tab6.4[, 5]=c(0, 4, 5, 1, 1, 0, 1, 0, 1, 0)
tab6.4[, 6]=c(0, 2, 3, 2, 0, 1, 0, 0, 0, 1)
tab6.4[, 7]=c(0, 0, 0, 0, 1, 0, 0, 1, 0, 0)
tab6.4[, 8]=c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
tab6.4[, 9]=c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
tab6.4[, 10]=c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
tab6.4[, 11]=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


# Zelterman (2004),pg.156, Tabel 6.9
#
# n	f_n	m_n	0	1	2	3	4	5	6	7+
#
# 1	267	12	255	12	0	0	0	0	0	0
# 2	285	48	239	44	2	0	0	0	0	0
# 3	202	80	143	41	15	3	0	0	0	0
# 4	110	54	69	30	9	2	0	0	0	0
# 5	104	103	43	34	15	9	3	0	0	0
# 6	50	67	15	18	8	5	3	0	1	0
# 7	21	38	4	4	7	4	2	0	0	0
# 8	12	28	1	2	4	3	1	1	0	0


tab6.9=matrix(0, nrow=8, ncol=12)
colnames(tab6.9)=c("n", "f_n", "m_n", as.character(0:8))
tab6.9[, 1]=c(1, 2, 3, 4, 5, 6, 7, 8)
tab6.9[, 2]=c(267, 285, 202, 110, 104, 50, 21, 12)
tab6.9[, 3]=c(12, 48, 80, 54, 103, 67, 38, 28)
tab6.9[, 4]=c(255, 239, 143, 69, 43, 15, 4, 1)
tab6.9[, 5]=c(12, 44, 41, 30, 34, 18, 4, 2)
tab6.9[, 6]=c(0, 2, 15, 9, 15, 8, 7, 4)
tab6.9[, 7]=c(0, 0, 3, 2, 9, 5, 4, 3)
tab6.9[, 8]=c(0, 0, 0, 0, 3, 3, 2, 1)
tab6.9[, 9]=c(0, 0, 0, 0, 0, 0, 0, 1)
tab6.9[, 10]=c(0, 0, 0, 0, 0, 1, 0, 0)
tab6.9[, 11]=c(0, 0, 0, 0, 0, 0, 0, 0)



# Zelterman (2004), pg.174, Tabel 7.1
#
#n	f_n	m_n	0	1	2	3	4	5	6
#
# 1	48	12	36	12	0	0	0	0	0
# 2	23	9	15	7	1	0	0	0	0
# 3	17	19	5	7	3	2	0	0	0
# 4	7	5	3	3	1	0	0	0	0
# 6	5	15	1	0	1	1	1	0	1


tab7.1=matrix(0, nrow=5, ncol=10)
colnames(tab7.1)=c("n", "f_n", "m_n", as.character(0:6))
tab7.1[, 1]=c(1, 2, 3, 4, 6)
tab7.1[, 2]=c(48, 23, 17, 7, 5)
tab7.1[, 3]=c(12, 9, 19, 5, 15)
tab7.1[, 4]=c(36, 15, 5, 3, 1)
tab7.1[, 5]=c(12, 7, 7, 3, 0)
tab7.1[, 6]=c(0, 1, 3, 1, 1)
tab7.1[, 7]=c(0, 0, 2, 0, 1)
tab7.1[, 8]=c(0, 0, 0, 0, 1)
tab7.1[, 9]=c(0, 0, 0, 0, 0)
tab7.1[, 10]=c(0, 0, 0, 0, 1)



# Zelterman (2004), pg.219, Tabel 8.2
#
# n	f_n	m_n	0	1	2	3	4	5	6	7+
#
# 1	267	12	255	12	0	0	0	0	0	0
# 2	285	48	239	44	2	0	0	0	0	0
# 3	202	80	143	41	15	3	0	0	0	0
# 4	110	54	69	30	9	2	0	0	0	0
# 5	104	103	43	34	15	9	3	0	0	0
# 6	50	67	15	18	8	5	3	0	1	0
# 7	21	38	4	4	7	4	2	0	0	0
# 8	12	28	1	2	4	3	1	1	0	0

tab8.2=matrix(0, nrow=8, ncol=12)
colnames(tab8.2)=c("n", "f_n", "m_n", as.character(0:8))
tab8.2[, 1]=c(1, 2, 3, 4, 5, 6, 7, 8)
tab8.2[, 2]=c(267, 285, 202, 110, 104, 50, 21, 12)
tab8.2[, 3]=c(12, 48, 80, 54, 103, 67, 38, 28)
tab8.2[, 4]=c(255, 239, 143, 69, 43, 15, 4, 1)
tab8.2[, 5]=c(12, 44, 41, 30, 34, 18, 4, 2)
tab8.2[, 6]=c(0, 2, 15, 9, 15, 8, 7, 4)
tab8.2[, 7]=c(0, 0, 3, 2, 9, 5, 4, 3)
tab8.2[, 8]=c(0, 0, 0, 0, 3, 3, 2, 1)
tab8.2[, 9]=c(0, 0, 0, 0, 0, 0, 0, 1)
tab8.2[, 10]=c(0, 0, 0, 0, 0, 1, 0, 0)
tab8.2[, 11]=c(0, 0, 0, 0, 0, 0, 0, 0)


# Zelterman (2004),pg. 245, Tabel 9.1.1
#
# n	f_n	m_n	0	1	2	3	4	5	6	7
#
# 5	71	60	30	27	9	5	0	0	0	0
# 6	156	95	86	51	14	4	1	0	0	0
# 7	224	163	111	73	31	8	1	0	0	0
# 8	150	104	79	44	23	3	0	1	0	0
# 9	70	48	32	29	8	1	0	0	0	0
# 10	12	9	5	5	2	0	0	0	0	0	

tab9.1.1=matrix(0, nrow=6, ncol=14)
colnames(tab9.1.1)=c("n", "f_n", "m_n", as.character(0:10))
tab9.1.1[, 1]=c(5, 6, 7, 8, 9, 10)
tab9.1.1[, 2]=c(71, 156, 224, 150, 70, 12)
tab9.1.1[, 3]=c(60, 95, 163, 104, 48, 9)
tab9.1.1[, 4]=c(30, 86, 111, 79, 32, 5)
tab9.1.1[, 5]=c(27, 51, 73, 44, 29, 5)
tab9.1.1[, 6]=c(9, 14, 31, 23, 8, 2)
tab9.1.1[, 7]=c(5, 4, 8, 3, 1, 0)
tab9.1.1[, 8]=c(0, 1, 1, 0, 0, 0)
tab9.1.1[, 9]=c(0, 0, 0, 1, 0, 0)
tab9.1.1[, 10]=c(0, 0, 0, 0, 0, 0)
tab9.1.1[, 11]=c(0, 0, 0, 0, 0, 0)


# Zelterman (2004),pg. 245, Tabel 9.1.2
#
# n	f_n	m_n	0	1	2	3	4	5	6	7
#
# 5	121	172	27	41	32	17	4	0	0	0
# 6	170	284	28	47	59	28	6	1	1	0
# 7	186	310	31	61	54	20	19	1	0	0
# 8	99	183	12	32	24	22	8	1	0	0
# 9	24	51	1	6	9	6	1	1	0	0
# 10	4	4	1	2	1	0	0	0	0	0


tab9.1.2=matrix(0, nrow=6, ncol=14)
colnames(tab9.1.2)=c("n", "f_n", "m_n", as.character(0:10))
tab9.1.2[, 1]=c(5, 6, 7, 8, 9, 10)
tab9.1.2[, 2]=c(121, 170, 186, 99, 24, 4)
tab9.1.2[, 3]=c(172, 284, 310, 183, 51, 4)
tab9.1.2[, 4]=c(27, 28, 31, 12, 1, 1)
tab9.1.2[, 5]=c(41, 47, 61, 32, 6, 2)
tab9.1.2[, 6]=c(32, 59, 54, 24, 9, 1)
tab9.1.2[, 7]=c(17, 28, 20, 22, 6, 0)
tab9.1.2[, 8]=c(4, 6, 19, 8, 1, 0)
tab9.1.2[, 9]=c(0, 1, 1, 1, 1, 0)
tab9.1.2[, 10]=c(0, 1, 0, 0, 0, 0)
tab9.1.2[, 11]=c(0, 0, 0, 0, 0, 0)


# Zelterman (2004),pg. 245, Tabel 9.1.3
#
# n	f_n	m_n	0	1	2	3	4	5	6	7
#
# 5	160	335	16	32	48	49	15	0	0	0
# 6	153	361	7	35	45	37	20	9	0	0
# 7	120	322	5	22	27	36	17	9	3	1
# 8	45	142	1	4	12	11	8	7	0	2
# 9	7	24	0	0	2	2	2	0	1	0
# 10	1	7	0	0	0	0	0	0	0	1



tab9.1.3=matrix(0, nrow=6, ncol=14)
colnames(tab9.1.3)=c("n", "f_n", "m_n", as.character(0:10))
tab9.1.3[, 1]=c(5, 6, 7, 8, 9, 10)
tab9.1.3[, 2]=c(160, 153, 120, 45, 7, 1)
tab9.1.3[, 3]=c(335, 361, 322, 142, 24, 7)
tab9.1.3[, 4]=c(16, 7, 5, 1, 0, 0)
tab9.1.3[, 5]=c(32, 35, 22, 4, 0, 0)
tab9.1.3[, 6]=c(48, 45, 27, 12, 2, 0)
tab9.1.3[, 7]=c(49, 37, 36, 11, 2, 0)
tab9.1.3[, 8]=c(15, 20, 17, 8, 2, 0)
tab9.1.3[, 9]=c(0, 9, 9, 7, 0, 0)
tab9.1.3[, 10]=c(0, 0, 3, 0, 1, 0)
tab9.1.3[, 11]=c(0, 0, 1, 2, 0, 1)



#Chung, C. S. and Brown, K. S. (1970).
#pg.632, Table2
#
# n	f_n	m_n	1	2	3	4	5+
#
# 1	73	73	73	0	0	0	0
# 2	120	120	89	31	0	0	0
# 3	105	105	84	20	1	0	0
# 4	56	56	38	17	1	0	0
# 5	31	31	20	8	3	0	0
# 6	18	18	12	2	4	0	0
# 7	6	6	5	0	0	1	0
# 8	11	11	8	1	1	1	0
# 9	8	8	6	2	0	0	0
# 10	2	2	1	1	0	0	0
# 11	1	1	1	0	0	0	0
# 12	1	1	0	0	0	1	0

tab2.cb=matrix(0, nrow=12, ncol=15)
colnames(tab2.cb)=c("n", "f_n", "m_n", as.character(1:5))
tab2.cb[, 1]=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
tab2.cb[, 2]=c(73, 120, 105, 56, 31, 18, 6, 11, 8, 2, 1, 1)
tab2.cb[, 3]=c(73, 120, 105, 56, 31, 18, 6, 11, 8, 2, 1, 1)
tab2.cb[, 4]=c(73, 89, 84, 38, 20, 12, 5, 8, 6, 1, 1, 0)
tab2.cb[, 5]=c(0, 31, 20, 17, 8, 2, 0, 1, 2, 1, 0, 0)
tab2.cb[, 6]=c(0, 0, 1, 1, 3, 4, 0, 1, 0, 0, 0, 0)
tab2.cb[, 7]=c(0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1)
tab2.cb[, 8]=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
print(tab2.cb)
tab2.cb[, 1]=tab2.cb[, 1]-1    # modýfy the sibling size by reducýng 1.
print(tab2.cb)




#Chung, C. S. and Brown, K. S. (1970).
#pg.633, Table3
#
# n	f_n	m_n	1	2	3	4	5+
#
# 1	8	8	8	0	0	0	0
# 2	14	14	12	2	0	0	0
# 3	16	16	12	4	0	0	0
# 4	13	13	7	5	0	1	0
# 5	20	20	15	4	1	0	0
# 6	11	11	4	5	1	1	0
# 7	7	7	3	2	2	0	0
# 8	7	7	2	1	3	1	0
# 9	8	8	5	2	0	1	0
# 10	3	3	2	0	0	1	0
# 11	3	3	1	2	0	0	0
# 13	1	1	1	0	0	0	0



tab3.cb=matrix(0, nrow=12, ncol=16)
colnames(tab3.cb)=c("n", "f_n", "m_n", as.character(1:5))
tab3.cb[, 1]=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13)
tab3.cb[, 2]=c(8, 14, 16, 13, 20, 11, 7, 7, 8, 3, 3, 1)
tab3.cb[, 3]=c(8, 14, 16, 13, 20, 11, 7, 7, 8, 3, 3, 1)
tab3.cb[, 4]=c(8, 12, 12, 7, 15, 4, 3, 2, 5, 2, 1, 1)
tab3.cb[, 5]=c(0, 2, 4, 5, 4, 5, 2, 1, 2, 0, 2, 0)
tab3.cb[, 6]=c(0, 0, 0, 0, 1, 1, 2, 3, 0, 0, 0, 0)
tab3.cb[, 7]=c(0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 0)
tab3.cb[, 8]=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
print(tab3.cb)
tab3.cb[, 1]=tab3.cb[, 1]-1    # modýfy the sibling size by reducýng 1.
print(tab3.cb)


#Chung, C. S. and Brown, K. S. (1970).
#pg.637, Table6
#
# n	f_n	m_n	0	1	2	3	4	5+
#
# 1	13	0	13	0	0	0	0	0
# 2	27	2	25	1	1	0	0	0
# 3	18	0	18	0	0	0	0	0
# 4	10	2	8	1	1	0	0	0
# 5	0	0	0	0	0	0	0	0
# 6	2	0	0	0	0	0	0	0
	
tab6.cb=matrix(0, nrow=6, ncol=10)
colnames(tab6.cb)=c("n", "f_n", "m_n", as.character(0:6))
tab6.cb[, 1]=c(1, 2, 3, 4, 5, 6)
tab6.cb[, 2]=c(13, 27, 18, 10, 0, 2)
tab6.cb[, 3]=c(0, 2, 0, 2, 0, 0)
tab6.cb[, 4]=c(13, 25, 18, 8, 0, 0)
tab6.cb[, 5]=c(0, 1, 0, 1, 0, 0)
tab6.cb[, 6]=c(0, 1, 0, 1, 0, 0)
tab6.cb[, 7]=c(0, 0, 0, 0, 0, 0)
tab6.cb[, 8]=c(0, 0, 0, 0, 0, 0)


#Chung, C. S. and Brown, K. S. (1970).
#pg.638, Table8
#
# n	f_n	m_n	0	1	2	3	4	5+
# 1	27	6	21	6	0	0	0	0
# 2	31	11	20	2	9	0	0	0
# 3	19	7	12	1	1	5	0	0
# 4	6	3	3	0	0	1	2	0
# 5	3	2	1	0	1	0	1	0
# 6	0	0	0	0	0	0	0	0
# 7	0	0	0	0	0	0	0	0
# 8	1	0	1	0	0	0	0	0


tab8.cb=matrix(0, nrow=8, ncol=12)
colnames(tab8.cb)=c("n", "f_n", "m_n", as.character(0:8))
tab8.cb[, 1]=c(1, 2, 3, 4, 5, 6, 7, 8)
tab8.cb[, 2]=c(27, 31, 19, 6, 3, 0, 0, 1)
tab8.cb[, 3]=c(6, 11, 7, 3, 2, 0, 0, 0)
tab8.cb[, 4]=c(21, 20, 12, 3, 1, 0, 0, 1)
tab8.cb[, 5]=c(6, 2, 1, 0, 0, 0, 0, 0)
tab8.cb[, 6]=c(0, 9, 1, 0, 1, 0, 0, 0)
tab8.cb[, 7]=c(0, 0, 5, 1, 0, 0, 0, 0)
tab8.cb[, 8]=c(0, 0, 0, 2, 1, 0, 0, 0)


# Pala et.al. (2008), Pg. 286.
#
# Table 2. Distributions and frequencies of Co-60 gamma-ray induced MN of pooled donors’ data.
#
#                         				Micronuclei distribution
# 
# Dose  Cells 	Total 	frequency 	   0   1    2     3     4     5     6    7
#			 MN
#	
#0    20.000 	259 		0.013 	19745 251 	4 	-	- 	- 	- 	-	
#0.1  20.000 	432 		0.021 	19579 410 	11 	-	- 	- 	- 	-
#0.25 20.000 	666 		0.033 	19354 628 	16 	2	- 	- 	- 	-
#0.5 	20.000 	1182 		0.059 	18851 1121 	23 	5 	- 	- 	- 	-
#0.75 20.000 	1660 		0.083 	18402 1549 	36 	13	- 	- 	- 	- 
#1 	20.000 	2376 		0.12 	17769 2119 	84	24 	3 	1 	- 	-
#2 	10.000 	2752 		0.28 	7449 	2373 	89 	48 	8 	5 	- 	- 
#3 	10.000 	4866 		0.49 	5812 	3773 	232	122	44 	15 	2 	- 
#4 	10.000 	7781 		0.78 	3663 	5481 	478 	222 	111 	37 	7 	1 
#5 	10.000 	10950 	1.10 		1652 	6818 	847 	420 	167 	70 	22 	4 

  


tab.pala1=matrix(0, nrow=10, ncol=13)
colnames(tab.pala1)=c("n", "f_n", "m_n", as.character(0:9))

tab.pala1[, 1]=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
tab.pala1[, 4]=c(19745, 19579, 19354, 18851, 18402, 17769, 7449, 5812, 3663, 1652)
tab.pala1[, 5]=c(251, 410, 628, 1121, 1549, 2119, 2373, 3773, 5481, 6818)
tab.pala1[, 6]=c(4, 11, 16, 23, 36, 84, 89, 232, 478, 847)
tab.pala1[, 7]=c(0, 0, 2, 5, 13, 24, 48, 122, 222, 420)
tab.pala1[, 8]=c(0, 0, 0, 0, 0, 3, 8, 44, 111, 167)
tab.pala1[, 9]=c(0, 0, 0, 0, 0, 1, 5, 15, 37, 70)
tab.pala1[, 10]=c(0, 0, 0, 0, 0, 0, 0, 2, 7, 22)
tab.pala1[, 11]=c(0, 0, 0, 0, 0, 0, 0, 0, 1, 4)


#  Elwood M. and Coldman A. (1981) Pg. 734
# Table 1. Diagnosis of breast cancer by number and sex of children
# 
# No of Children 	   0 	  1 	    2 	       3 	           4
#                 -----  -----  ----------   ----------    --------------
# No of boys 	   0 	 0  1    0   1  2    0  1  2  3     0  1  2  3  4
#
# No of patients    284 93 71   65 134 83   26 71 75 26    11 21 30 28  4


#                       number of boys

# n	f_n	m_n	0	1	2	3	4
# 0	289	0     289 
# 1	164	71	93	71
# 2	282   300	65	134	83
# 3	198	299	26	71	75	26
# 4	 94	181	11	21	30	28	4

tab.ec=matrix(0, nrow=5, ncol=8)
colnames(tab.ec)=c("n", "f_n", "m_n", as.character(0:4))
tab.ec[,1]=c(0,1,2,3,4)
tab.ec[,2]=c(289, 164, 282, 198, 94)
tab.ec[,3]=c(0, 71, 300, 299, 181)
tab.ec[,4]=c(289, 93, 65, 26, 11)
tab.ec[,5]=c(0, 71, 134, 71, 21)
tab.ec[,6]=c(0, 0, 83, 75, 30)
tab.ec[,7]=c(0, 0, 0, 26, 28)
tab.ec[,8]=c(0, 0, 0, 0, 4)

######################
# Rao C. R. (1997) pg. 115 
# Table 4.6 Distribution of birth rank s and family size n
#--------------------------------------------------
#   n= 1  2  3  4
# s
# 1  21  22  17  11
# 2  -   10  14  10    
# 3  -   -   9   13
# 4  -   -   -   13

tab4.6rao=matrix(0, nrow=4, ncol=7)
colnames(tab4.6rao)=c("n", "f_n", "m_n", as.character(1:4))
tab4.6rao[, 1]=c(1, 2, 3, 4)
tab4.6rao[, 4]=c(21, 22, 17, 11)
tab4.6rao[, 5]=c(0, 10, 14, 10)
tab4.6rao[, 6]=c(0, 0, 9, 13)
tab4.6rao[, 7]=c(0, 0, 0, 13)
tab4.6rao[, 1]=tab4.6rao[, 1]-1    # modify the sibling size by reducýng 1.

# Rao C. R. (1997) pg. 115 
# Table 4.7 Distribution of birth ranks 
#           and family size n<=4 
#           among staff members. 
#           (University of Pittsburgh)
#--------------------------------------------------
#   n= 1  2  3  4
# s
# 1  7  14  9  6
# 2  -  6   4  2    
# 3  -   -  2  0
# 4  -   -  -  0

tab4.7rao=matrix(0, nrow=4, ncol=7)
colnames(tab4.7rao)=c("n", "f_n", "m_n", as.character(1:4))
tab4.7rao[, 1]=c(1, 2, 3, 4)
tab4.7rao[, 4]=c(7, 14, 9, 11)
tab4.7rao[, 5]=c(0, 6, 4, 2)
tab4.7rao[, 6]=c(0, 0, 2, 0)
tab4.7rao[, 7]=c(0, 0, 0, 0)
tab4.7rao[, 1]=tab4.7rao[, 1]-1    # modýfy the sibling size by reducýng 1.

# Rao C. R. (1997) pg. 155 
# Table 5.7 Mean scores on the National Merit
#           Scholarship Qualification Test,1965, 
#           by place in family configuration in USA
#
# Family 				Birth order
# size	1	    2	     3	 4	  5	
# 1	    103.76      -	     -	 -	  -
# 2	    106.21   104.44    -       -      -
# 3       106.14   102.89  102.71    -      -
# 4	    105.59   103.05  101.30   100.18  -
# 5	    104.39   101.71  99.37    97.69   96.87

tab5.7rao=matrix(0, nrow=5, ncol=8)
colnames(tab5.7rao)=c("n", "f_n", "m_n", as.character(1:5))
tab5.7rao[, 1]=c(1, 2, 3, 4, 5)
tab5.7rao[, 4]=c(103.76, 106.21, 106.14, 105.59, 104.39)
tab5.7rao[, 5]=c(0, 104.44, 102.89, 103.05, 101.71)
tab5.7rao[, 6]=c(0, 0, 102.71, 101.30, 99.37)
tab5.7rao[, 7]=c(0, 0, 0, 100.18, 97.69)
tab5.7rao[, 8]=c(0, 0, 0,  0, 96.87)
tab5.7rao[, 1]=tab5.7rao[, 1]-1    # modify the sibling size by reducýng 1.



# Hubertlacey J., Gowers S. G. and Bhat A. V. (1991) pp.492
# Table 1. Distribution of birth order
#
# Sibship  Total no.     Patients order within sibship
# size    of families    1 	2 	3 	4 	5
#
#  1		14		14		
#  2		81		49	32
#  3		62		20	14	28
#  4		44		10	12	14	8
#  5		24		10	 2	 0	4	8
# Total	225		103	60	42	12	8


tab.hgb=matrix(0, nrow=5, ncol=8)
colnames(tab.hgb)=c("n", "f_n", "m_n", as.character(1:5))
tab.hgb[, 1]=c(1, 2, 3, 4, 5)
tab.hgb[, 2]=c(14, 81, 62, 44, 24)
tab.hgb[, 4]=c(14, 49, 20, 10, 10)
tab.hgb[, 5]=c(0, 32, 14, 12, 2)
tab.hgb[, 6]=c(0, 0, 28, 14, 0)
tab.hgb[, 7]=c(0, 0, 0, 8, 4)
tab.hgb[, 8]=c(0, 0, 0, 0, 8)
tab.hgb[, 1]=tab.hgb[, 1]-1    # modify the sibling size by reducing 1.



# Haukka, J. K., Suvisaari J. and Lönnqvist J. (2004) (second page of paper)
#
# Table 3: A family study of schizophrenia in Finland 
# in birth cohorts born from 1950 to 1976. 
# Family structure of families (N = 2605)
# with at least one sibling with schizophrenia; 
# also siblings born before 1950 and after 1976 are included.
# 
# N. of cases		Family size
# in  the 			
# family 	
#		 2	3	4	5	6	7	8	9	10	11	12
#
# 1		3289	2012	984	369	171	62	22	9	5	1	0
# 2 		282	271	190	82	35	16	12	5	1	0	0
# 3		0	22	26	11	10	8	4	0	1	0	1
# 4		0	0	4	2	1	1	1	0	0	0	0
# 5		0	0	0	1	0	1	0	0	0	0	0
# 6		0	0	0	0	0	0	1	1	1	0	0
#
#






tab.hauk1=matrix(0, nrow=6, ncol=14)
colnames(tab.hauk1)=c("n", "f_n", "m_n", as.character(2:12))
tab.hauk1[, 1]=c(1, 2, 3, 4, 5, 6)
tab.hauk1[, 4]=c(3289, 282, 0, 0, 0, 0)
tab.hauk1[, 5]=c(2012, 271, 22, 0, 0, 0)
tab.hauk1[, 6]=c(984, 190, 26, 4, 0, 0)
tab.hauk1[, 7]=c(369, 82, 11, 2, 1, 0)
tab.hauk1[, 8]=c(171, 35, 10, 1, 0, 0)
tab.hauk1[, 9]=c(62, 16, 8, 1, 1, 0)
tab.hauk1[, 10]=c(22, 12, 4, 1, 0, 1)
tab.hauk1[, 11]=c(9, 5, 0, 0, 0, 1)
tab.hauk1[, 12]=c(5, 1, 1, 0, 0, 1)
tab.hauk1[, 13]=c(1, 0, 0, 0, 0, 0)
tab.hauk1[, 14]=c(0, 0, 1, 0, 0, 0)
tab.hauk1[, 1]=tab.hauk1[, 1]-1    # modify the sibling size by reducing 1.

# col as rows

tab.hauk2=matrix(0, nrow=11, ncol=15)
colnames(tab.hauk2)=c("n", "f_n", "m_n", as.character(1:6))
tab.hauk2[, 1]=c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
tab.hauk2[, 4]=c(3289, 2012, 984, 369, 171, 62, 22, 9, 5, 1, 0)
tab.hauk2[, 5]=c(282, 271, 190, 82, 35, 16, 12, 5, 1, 0, 0)
tab.hauk2[, 6]=c(0, 22, 26, 11, 10, 8, 4, 0, 1, 0, 1)
tab.hauk2[, 7]=c(0, 0, 4, 2, 1, 1, 1, 0, 0, 0, 0)
tab.hauk2[, 8]=c(0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0)
tab.hauk2[, 9]=c(0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0)
tab.hauk2[, 1]=tab.hauk2[, 1]-2    # modify the sibling size by reducing 2.


# Omariba D. W. R., Rajulton F., Beaujot R. (2008) pp.320.
# Table 1: Distribution of children and infant deaths in Kenya, DHS1998
#						Infant deaths in the family
#
# Children in famuly	0	1	2	3	4	5
#
# 1			  1038	62	0	0	0	0
# 2			   897	80	4	0	0	0
# 3			   666	105	9	0	0	0
# 4			   575	97	18	3	1	0
# 5			   406	92	21	8	1	0
# 6			   372	105	27	6	3	1
# 7			   232	90	17	16	1	0
# 8			   157	72	25	12	5	2
# 9			   131	57	20	9	2	3
# 10-15		   92		68	42	22	17	6


tab.omariba=matrix(0, nrow=10, ncol=14)
colnames(tab.omariba)=c("n", "f_n", "m_n", as.character(0:5))
tab.omariba[,1]=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
tab.omariba[,4]=c(1038, 897, 666, 575, 406, 372, 232, 157, 131, 92)
tab.omariba[,5]=c(62, 80, 105, 97, 92, 105, 90, 72, 57, 68)
tab.omariba[,6]=c(0, 4, 9, 18, 21, 27, 17, 25, 20, 42)
tab.omariba[,7]=c(0, 0, 0, 3, 8, 6, 16, 12, 9, 22)
tab.omariba[,8]=c(0, 0, 0, 1, 1, 3, 1, 5, 2, 17)
tab.omariba[,9]=c(0, 0, 0, 0, 0, 1, 0, 2, 3, 6)



# Resource: http://www.statsci.org/data/general/fetaimpl.html
#
# Faddy, M. J., and Bosch, R. J. (1999). pp. 621
#
#Table 1.Frequency of observed numbers of fetal implants of mice
#exposed to varying doses of 2, 4, 5-T.

# In an experiment where pregnant mice were exposed to 
# the herbicide 2,4,5-T (the active component in Agent Orange), 
# the number of fetal implants in utero were recorded. 
# The data give the frequency distribution of implants at each of 
# seven dose levels measured in mg/kg of body weight.
#
# On days 6-14 after mating, pregnant dams were dosed by 
# gavage with one of the doses of 2,4,5-T. Prior to giving birth, 
# the dams were sacrificed and the number of viable, dead 
# and reabsorbed foetuses in the uterus of the dam were determined. 
# The data here gives the number of surviving viable implants. 
# An outcome of zero implants cannot be distinguished from
# a non-pregnant outcome so any zero implant outcomes were excluded.
#
# 
#--------------------------------------------------------------------------------
# 
# Variable  Descriptions 
# 
#--------------------------------------------------------------------------------
# 
# Dose  Dose of 2,4,5-T in mg/kg/day 
# Implants  Number of surviving implants 
# Frequency  Number of mice with that number of implants 
# 
# --------------------------------------------------------------------------------
# Comment: The counts are underdispersed and left skew relative to the Poisson distribution.
#
# This data is also used by Holson et. al (1992)
# Holson, J.F., Gaines, T.B., Nelson, C.J., LaBorde, J.B., Gaylor, D.W., Sheehan, D.M. and Young, J.F. (1992). 
# Developmental toxicity of 2,4,5-trichlorophenoxyacetic acid (2,4,5-T):
# I. Multireplicated dose-response studies in four inbred strains and one outbred stock of mice. 
# Fundamental and Applied Toxicology 19, 286-297.

tab.fad=matrix(0, nrow=7, ncol=21)

tab.fad[1,]=c(6, 11, 4, 9, 9, 11, 17, 28, 41, 69, 83, 120, 112, 80, 52, 25, 15, 6, 0, 0, 0)
tab.fad[2,]=c(5, 2, 4, 8, 5, 6, 10, 9, 10, 14, 30, 45, 56, 38, 31, 27, 4, 1, 2, 0, 0)
tab.fad[3,]=c(6, 7, 9, 11, 10, 15, 19, 22, 36, 78, 101, 99, 133, 82, 52, 27, 9, 3, 1, 0, 2) 
tab.fad[4,]=c(1, 0, 2, 1, 2, 2, 2, 1, 8, 4, 20, 20, 13, 12, 7, 1, 1, 1, 0, 0, 0)
tab.fad[5,]=c(4, 3, 6, 8, 10, 13, 17, 20, 31, 48, 85, 91, 111, 67, 40, 20, 13, 5, 0, 0, 0)
tab.fad[6,]=c(0, 0, 0, 0, 2, 0, 0, 3, 7, 4, 8, 8, 8, 3, 1, 0, 0, 0, 0, 0, 0)
tab.fad[7,]=c(0, 1, 3, 1, 1, 1, 2, 6, 3, 11, 9, 15, 14, 11, 3, 2, 0, 0, 0, 0, 0)


# Ekholm A., Smith P. W. F. and McDonald ()
# Marginal regression analysis of a multivariate binary response
# Biometrika (1995), 82,4, pp. 847-54
# Table 1. Observed and fitted frequencies of
#  impaired pulmonary function 
#  among 203 siblings from 100 families
#
# Comment:I realize that this data is the Table 6.1 in Zelterman.

tab1.ekholm=matrix(0, nrow=6, ncol=10)
colnames(tab1.ekholm)=c("n", "f_n", "m_n", as.character(0:6))
tab1.ekholm[, 1]=c(1, 2, 3, 4, 5, 6)
tab1.ekholm[, 2]=c(48, 23, 17, 7, 0, 5)
tab1.ekholm[, 3]=c(12, 9, 19, 5, 0, 15)
tab1.ekholm[, 4]=c(36, 15, 5, 3, 0, 1)
tab1.ekholm[, 5]=c(12, 7, 7, 3, 0, 0)
tab1.ekholm[, 6]=c(0, 1, 3, 1, 0, 1)
tab1.ekholm[, 7]=c(0, 0, 2, 0, 0, 1)
tab1.ekholm[, 8]=c(0, 0, 0, 0, 0, 1)
tab1.ekholm[, 9]=c(0, 0, 0, 0, 0, 0)
tab1.ekholm[, 10]=c(0, 0, 0, 0, 0, 1)


# Comment:I realize that this data is the Table 6.1 in Zelterman.





# Calculate dispersion for each row and accumulated rows.

source("analysis_WordCount.r")
source("functions-altham-qwbinom.r")


rd=function(x=freq5)
{
# Compute sample dispersion index relative to the binomial 
# distribution of the thea same mean and support.
#
  M=length(x)-1
  Support=0:M
  n=sum(x)               # sample size
  Mean=sum(x*Support)/n
  Var=(sum(x*Support^2)-n*Mean^2)/(n-1)
  
  p.hat=Mean/M
  rd=Var/(Mean*(1-p.hat))

  return(rd)
}



rd.by.row=function(x, row0=3)
{
# Calculate rd for rows bigger or equal to row0.
#  
  rd0=NULL
  m=nrow(x)
  for(i in row0:m)
  {
    y=x[i, 4:(4+x[i, 1])]
#    print(y)
#   readline()
    rd0=c(rd0, rd(y))
  }
  print(rd0)

# Calculate rd for accumulated rows bigger or equal to row0.
#
  rd1=NULL
  for(i in row0:m)
  {
    tmp=x[1:i, 4:(4+x[i, 1])]
    y=apply(tmp, 2, sum)
#    print(y)
#    readline()
    rd1=c(rd1, rd(y))
  }  
  print(rd1)

}




x=tab6.1; row0=2
rd.by.row(x, row0)
#[1] 1.405016 0.973913 3.333333
#[1] 1.246305 1.163826 1.795735

x=tab6.2; row0=3
rd.by.row(x, row0)
#[1] 1.0666667 1.2571429 0.0000000       NaN       NaN       NaN 0.6666667
#[8]       NaN
#[1] 1.3933333 1.1884754 0.9129968 0.8016807 0.7264038 0.6698565 0.7007992
#[8] 0.7450817

x=tab6.3; row0=3
rd.by.row(x, row0)
#[1] 1.393973 1.224720 1.479030 1.749551 1.164075 1.136746
#[1] 1.265661 1.270871 1.464611 1.602024 1.640585 1.688020

x=tab6.4; row0=2
rd.by.row(x, row0)
#[1] 1.142857 1.055195 1.333333 4.079861 2.352941      NaN      NaN      NaN
#[9]      NaN
#[1] 1.2375000 0.9748840 0.9271978 1.5357025 2.2633528 2.1255080 2.0181461
#[8] 1.9347440 1.7738370


x=tab6.9; row0=2
rd.by.row(x, row0)
#[1] 1.002560 1.393973 1.224720 1.479030 1.749551 1.164075 1.136746
#[1] 1.014865 1.265661 1.270871 1.464611 1.602024 1.640585 1.688020

x=tab7.1; row0=2
rd.by.row(x, row0)
#[1] 1.079989 1.405016 0.973913 3.333333
#[1] 0.9516163 1.2463054 1.1638263 1.7957351

x=tab8.2; row0=2
rd.by.row(x, row0)
#[1] 1.002560 1.393973 1.224720 1.479030 1.749551 1.164075 1.136746
#[1] 1.014865 1.265661 1.270871 1.464611 1.602024 1.640585 1.688020


x=tab9.1.1; row0=2
rd.by.row(x, row0)
#[1] 1.1926253 1.1443222 1.2282897 0.8484611 0.8190008
#[1] 1.178476 1.148752 1.154472 1.110946 1.095140

x=tab9.1.2; row0=2
rd.by.row(x, row0)
#[1] 1.0679015 1.1464865 1.0102857 0.7667054 0.7407407
#[1] 1.0824836 1.0761740 1.0361686 0.9998482 0.9777143

x=tab9.1.3; row0=2
rd.by.row(x, row0)

#[1] 1.1085549 1.2389633 1.2358891 0.9198718       NaN
#[1] 1.031345 1.056501 1.045275 1.000595 0.983320

x=tab2.cb; row0=3
#[1] 0.9939927 0.8793552 1.1380087 1.4823529 3.2727273 2.1328638 0.8847926
# [8] 1.0588235       NaN       NaN
# [1] 0.9469910 0.9167310 0.9481779 0.9992552 1.0449792 1.0964377 1.0854423
# [8] 1.0776435 1.0750877 1.1139362

x=tab3.cb; row0=4
rd.by.row(x, row0)
#[1] 1.5463710 1.1759128 1.1977778 1.1018519 1.1307692 1.9525424 3.3750000
#[8] 0.5357143       NaN
#[1] 1.295704 1.220216 1.275200 1.248296 1.319622 1.349624 1.401066 1.357442
#[9] 1.349597

x=tab6.cb; row0=4
rd.by.row(x, row0)
#[1] 1.641642      NaN      NaN
#[1] 1.638125 1.630768 1.625900

x=tab8.cb; row0=3
#[1] 2.823837 4.195804 3.333333      NaN      NaN      NaN
#[1] 1.958591 2.173878 2.159591 2.101686 2.062190 2.040023

x=tab.pala1; row0=3
rd.by.row(x, row0)

# [1] 1.0325715 1.0254494 1.0287586 1.0613517 1.0111863 1.0121295 0.9024853
# [8] 0.7700856
# [1] 1.035045 1.038338 1.040387 1.058412 1.083806 1.155894 1.241475 1.324796

x=tab.ec; row0=2
rd.by.row(x, row0)
# [1] 1.006135 1.049479 1.038827 1.168057
# [1] 1.002212 1.262624 1.298166 1.326066

x=tab4.6rao; row0=2
rd.by.row(x, row0)
# [1] 1.032258 1.303419 1.726430
# [1] 1.019231 1.275621 1.685135

x=tab4.7rao; row0=2
rd.by.row(x, row0)
# [1] 1.0526316 1.4123377 0.9662162
# [1] 1.038462 1.170732 1.081749

x=tab5.7rao; row0=2
rd.by.row(x, row0)

# [1] 1.004770 1.344130 1.675084 2.012060
# [1] 1.003191 1.252733 1.499520 1.745202

x=tab.hgb; row=2
rd.by.row(x, row0)

# [1] 1.012500 1.565999 1.455882 3.389943
# [1] 1.010638 1.340577 1.344202 1.637819

x=tab.hauk1; row=2
rd.by.row(x, row0)

#[1] 1.0018116 0.7176538 0.5142857       NaN       NaN
#[1] 1.000171 1.266282 1.426001 1.573328 1.658122


x=tab.hauk2; row=2
rd.by.row(x, row0)
# [1] 1.000438 1.127771 1.166519 1.271419 1.585573 1.758242 2.763158 3.591837
# [9]      NaN      NaN
# [1] 1.000171 1.058761 1.078809 1.092004 1.107626 1.126459 1.138160 1.152862
# [9] 1.150600 1.149903

x=tab.omariba; row=2
rd.by.row(x, row0)
#[1] 1.049289 1.044839 1.312501 1.411965 1.470317 1.416150 1.638758 1.679281
#[9] 1.645886
#[1] 1.018429 1.033024 1.135476 1.226816 1.307125 1.352185 1.439626 1.489120
#[9] 1.653125

#Comments:
#
#tab6.2 and tab3.cb for underdispersion overdispersion in the same data 


x=tab1.ekholm; row0=2
rd.by.row(x, row0)

#[1] 1.079989 1.405016 0.973913      NaN 3.333333
#[1] 0.9516163 1.2463054 1.1638263 1.1333773 1.7957351

x=tab.fad;
ds=matrix(NA,nrow=1, ncol=7)

for(i in 1:7)
{
ds[,i]=rd(x[i,]) # dispersion of each row
}
ds

#> ds
#         [,1]     [,2]     [,3]     [,4]     [,5]      [,6]     [,7]
#[1,] 1.989590 2.513324 1.995054 1.776459 1.870784 0.9949466 1.853231








