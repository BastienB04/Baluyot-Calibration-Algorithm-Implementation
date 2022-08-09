//Implementation of Baluyot Calibration Algorithm
    //written by Bastien Carreon Baluyot (GitHub: BastienB04)

//This is a simple implementation of the algorithm on C++
//This implementation is mainly for testing the algorithm's correctness.

//This implementation takes in 2 inputs and uses the original context of "SoC" as the input, though not necessarily capped from 0-100%
//This implementation has 3 options for the input: Manual, Automated (with Pseudorandom Number Generation) or Automated (by reading values from external .txt files)
//All necessary data is logged and output to external .txt files when the program terminates.


#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

//Function to read input .txt file line by line:
bool getFileContent(std::string fileName, std::vector<double>& v_dbl){
    std::ifstream in(fileName.c_str());

    if(!in){
      std::cerr << "Cannot open the File : "<<fileName<<std::endl;
      return false;
    }
  
    std::string str;

    while (std::getline(in, str)){
        if(str.size() > 0){
          v_dbl.push_back(std::stoi(str));
        }
    }

    in.close();
    return true;
}

//Function to calculate MSE of parameter SoC-t dataset from 'true' dataset (SoC(t) = 10t):
double MSE(std::vector<double>& v, std::vector<double>& time){
  double sum = 0, negatives = 0;
  
  for (int i=0; i<time.size(); i++){
    if (!(v[i] < 0)){  //Disregard negative elements, which could be present in 'unclean' calibrated dataset
      sum = sum + pow(std::abs(v[i] - (10*i)),2);
    }
    else{
      negatives++;
    }
  }

  return sum/(v.size() - negatives);
}

//Function to remove invalid SoC measurements from calibrated dataset:
void cleanVector (std::vector<double>& SoC_calibrated, std::vector<double>& time_vector){

  std::vector<double> temp1, temp2;
  
  for (int i=0; i<SoC_calibrated.size(); i++){
    if (!(SoC_calibrated[i] < 0)){
      temp1.push_back(SoC_calibrated[i]);
      temp2.push_back(time_vector[i]);    
    }
  }
  
  SoC_calibrated = temp1;
  time_vector = temp2;
}


//Function to log all relevant data to external .txt files:
void outputData (std::vector<double>& time_vector, std::vector<double>& SoC1, std::vector<double>& SoC1_fit, std::vector<double>& SoC2, std::vector<double>& SoC2_fit, std::vector<double>& Delta, std::vector<double>& Delta_fit, std::vector<double>& SoC_calibrated, std::vector<double>& EPSILON_vector, std::vector<double>& MSE_SoC1, std::vector<double>& MSE_SoC2, std::vector<double>& MSE_SoC_calibrated){
  
      std::ofstream file_out1;
      file_out1.open("Time.txt");
      for (int i=0; i<time_vector.size(); i++){
        file_out1 << time_vector[i] << std::endl;
      }
      file_out1.close();
  
      
      std::ofstream file_out2;
      file_out2.open("SoC1.txt");
      for (int i=0; i<SoC1.size(); i++){
        file_out2 << SoC1[i] << std::endl;
      }
      file_out2.close();


      std::ofstream file_out3;
      file_out3.open("SoC1 Fit.txt");
      for (int i=0; i<SoC1_fit.size(); i++){
        file_out3 << SoC1_fit[i] << std::endl;
      }
      file_out3.close();


      std::ofstream file_out4;
      file_out4.open("SoC2.txt");
      for (int i=0; i<SoC2.size(); i++){
        file_out4 << SoC2[i] << std::endl;
      }
      file_out4.close();


      std::ofstream file_out5;
      file_out5.open("SoC2 Fit.txt");
      for (int i=0; i<SoC2_fit.size(); i++){
        file_out5 << SoC2_fit[i] << std::endl;
      }
      file_out5.close();


      std::ofstream file_out6;
      file_out6.open("Delta.txt");
      for (int i=0; i<Delta.size(); i++){
        file_out6 << Delta[i] << std::endl;
      }
      file_out6.close();


      std::ofstream file_out7;
      file_out7.open("Delta fit.txt");
      for (int i=0; i<Delta_fit.size(); i++){
        file_out7 << Delta_fit[i] << std::endl;
      }
      file_out7.close();


      std::ofstream file_out8;
      file_out8.open("SoC Calibrated.txt");
      for (int i=0; i<SoC_calibrated.size(); i++){
        file_out8 << SoC_calibrated[i] << std::endl;
      }
      file_out8.close();



      
      std::ofstream file_out9;
      file_out9.open("Epsilon.txt");
      for (int i=0; i<EPSILON_vector.size(); i++){
        file_out9 << EPSILON_vector[i] << std::endl;
      }
      file_out9.close();


  
      std::ofstream file_out10;
      file_out10.open("MSE SoC1.txt");
      for (int i=0; i<MSE_SoC1.size(); i++){
        file_out10 << MSE_SoC1[i] << std::endl;
      }
      file_out10.close();

      std::ofstream file_out11;
      file_out11.open("MSE SoC2.txt");
      for (int i=0; i<MSE_SoC2.size(); i++){
        file_out11 << MSE_SoC2[i] << std::endl;
      }
      file_out11.close();

      std::ofstream file_out12;
      file_out12.open("MSE SoC Calibrated.txt");
      for (int i=0; i<MSE_SoC_calibrated.size(); i++){
        file_out12 << MSE_SoC_calibrated[i] << std::endl;
      }
      file_out12.close();

      cleanVector(SoC_calibrated,time_vector);
  
      std::ofstream file_out13;
      file_out13.open("SoC Calibrated (CLEANED).txt");
      for (int i=0; i<SoC_calibrated.size(); i++){
        file_out13 << SoC_calibrated[i] << std::endl;
      }
      file_out13.close();

      std::ofstream file_out14;
      file_out14.open("Time Calibrated (CLEANED).txt");
      for (int i=0; i<time_vector.size(); i++){
        file_out14 << time_vector[i] << std::endl;
      }
      file_out14.close();
}





int main() {

  std::cout << "START PROGRAM" << std::endl;
  
  //GLOBAL VARIABLES:
  
  double EPSILON = 1, EPSILON_step = 0.0005;
  double time_step = 1;
  bool Charging = true, Discharging = false; //this can be changed for testing purposes but should be preset depending on the purpose.


  std::vector<double> time_vector, SoC1, SoC2, SoC_calibrated, MSE_SoC1, MSE_SoC2, MSE_SoC_calibrated, EPSILON_vector;
  time_vector.push_back(0);
  double next_SoC1, next_SoC2, next_SoC_calibrated;

    //Curve fitting variables:
  double xsum_SoC1, x2sum_SoC1, ysum_SoC1, xysum_SoC1,
  xsum_SoC2, x2sum_SoC2, ysum_SoC2, xysum_SoC2,
  xsum_D, x2sum_D, ysum_D, xysum_D;
  double m1, c1, m2,c2, m_D, c_D;
  int n;

    //First validity criterion variables:
  double Delta_max, Delta_min, Delta_range;
  double SoC1_predicted, SoC2_predicted;
  double PvA_SoC1, PvA_SoC2;

    //Calibration variables:
  bool SoC1_INVALID = false, SoC2_INVALID = false;
  int ERRORMSG = 0;
  int weight1 = 1, weight2 = 1;
  
  int InvalidLoopCounter = 0, WhileLoopCounter = 0;

  



  //If automating SoC inputs via PRNG:
  int rand_min = 0;
  srand(time(NULL));



  //If automating SoC inputs via .txt file reading:
  std::vector<double> SoC1_auto, SoC2_auto;
  bool temp1, temp2;
  temp1 = getFileContent("SoC1_input.txt",SoC1_auto);
  temp2 = getFileContent("SoC2_input.txt",SoC2_auto);


  


  while (true){

    /*
    //Manual input:
    std::cout << "Next SoC1: " << std::endl;
    std::cin >> next_SoC1;
    std::cout << "Next SoC2: " << std::endl;
    std::cin >> next_SoC2;
    */

    //Automated input: (Randomly generated)
    next_SoC1 = rand() % 30 + rand_min;  //Change number to control degree of variance
    next_SoC2 = rand() % 30 + rand_min;  //Change number to control degree of variance
    rand_min = rand_min + 10;            //Change number to control gradient of 'true' SoC(t) function

    /*
    //Automated input: (From .txt file)
    next_SoC1 = SoC1_auto[WhileLoopCounter];
    next_SoC2 = SoC2_auto[WhileLoopCounter];
    */
    
    //First 3 measurements shouldn't trigger the calibration algorithm just yet, so continue:
    if (WhileLoopCounter <= 3){
      SoC1.push_back(next_SoC1);
      SoC2.push_back(next_SoC2);
      SoC_calibrated.push_back(-20);

      time_vector.push_back(time_vector.back()+time_step);
      EPSILON_vector.push_back(EPSILON);
      WhileLoopCounter++;

      continue;
    }


    //Create Delta vector from scratch each loop using the latest SoC1 and SoC2 vectors:
    std::vector<double> Delta;
    for (int i = 0; i < std::ceil((SoC1.size()+SoC2.size())/2); i++) {
      Delta.push_back(std::abs(SoC1[i] - SoC2[i]));
    }


   
  
    //Curve fit all SoC1 measurements in dataset using Linear Least Squares
    //(This would theoretically be modelling the SoC_n measurements for n = 1 and t = 0 to k)
    
    n = SoC1.size();
    xsum_SoC1 = 0;
    x2sum_SoC1 = 0;
    ysum_SoC1 = 0;
    xysum_SoC1 = 0;

    for (int i = 0; i < n; i++) {
      xsum_SoC1 = xsum_SoC1 + time_vector[i];
      ysum_SoC1 = ysum_SoC1 + SoC1[i];
      x2sum_SoC1 = x2sum_SoC1 + pow(time_vector[i], 2);
      xysum_SoC1 = xysum_SoC1 + SoC1[i] * time_vector[i];
    }

    m1 = (n * xysum_SoC1 - xsum_SoC1 * ysum_SoC1) / (n * x2sum_SoC1 - xsum_SoC1 * xsum_SoC1);
    c1 = (x2sum_SoC1 * ysum_SoC1 - xsum_SoC1 * xysum_SoC1) / (x2sum_SoC1 * n - xsum_SoC1 * xsum_SoC1);
  

    std::vector<double> SoC1_fit;
    for (int i = 0; i < n; i++) {
      SoC1_fit.push_back(m1 * time_vector[i] + c1);
    }
  

    std::cout << "The best fit line for SoC1(t) is: SoC1 = " << m1 << " t + " << c1 <<std::endl;
    std::cout << "The predicted next SoC1 according to the best fit line would have been: " << (m1 * (time_vector.back())) + c1 << std::endl;
    std::cout << "The predicted next SoC1 of the next time interval according to the best fit line would be: " << (m1 * (time_vector.back() + time_step)) + c1 << std::endl;
  
 
    //Curve fit all SoC2 measurements in dataset using Linear Least Squares
    //(This would theoretically be modelling the SoC_n measurements for n = 2 and t = 0 to k)
    
    n = SoC2.size();
    xsum_SoC2 = 0;
    x2sum_SoC2 = 0;
    ysum_SoC2 = 0;
    xysum_SoC2 = 0;
    
    for (int i = 0; i < n; i++) {
      xsum_SoC2 = xsum_SoC2 + time_vector[i];
      ysum_SoC2 = ysum_SoC2 + SoC2[i];
      x2sum_SoC2 = x2sum_SoC2 + pow(time_vector[i], 2);
      xysum_SoC2 = xysum_SoC2 + SoC2[i] * time_vector[i];
    }

    m2 = (n * xysum_SoC2 - xsum_SoC2 * ysum_SoC2) / (n * x2sum_SoC2 - xsum_SoC1 * xsum_SoC2);
    c2 = (x2sum_SoC2 * ysum_SoC2 - xsum_SoC2 * xysum_SoC2) / (x2sum_SoC2 * n - xsum_SoC2 * xsum_SoC2);
  

    std::vector<double> SoC2_fit;
    for (int i = 0; i < n; i++) {
      SoC2_fit.push_back(m2 * time_vector[i] + c2);
    }

    std::cout << "The best fit line for SoC2(t) is: SoC2 = " << m2 << " t + " << c2 <<std::endl;
    std::cout << "The predicted next SoC2 according to the best fit line would have been: " << (m2 * (time_vector.back())) + c2 << std::endl;
    std::cout << "The predicted next SoC2 of the next time interval according to the best fit line would be: " << (m2 * (time_vector.back() + time_step)) + c2 << std::endl;
  
  
  
  
    
    // Curve fit Delta
    n = Delta.size();
    xsum_D = 0;
    x2sum_D = 0;
    ysum_D = 0;
    xysum_D = 0;
    
    for (int i = 0; i < n; i++) {
      xsum_D = xsum_D + time_vector[i];
      ysum_D = ysum_D + Delta[i];
      x2sum_D = x2sum_D + pow(time_vector[i], 2);
      xysum_D = xysum_D + Delta[i] * time_vector[i];
    }
    
    m_D = (n * xysum_D - xsum_D * ysum_D) / (n * x2sum_D - xsum_D * xsum_D);

    c_D = (x2sum_D * ysum_D - xsum_D * xysum_D) / (x2sum_D * n - xsum_D * xsum_D);

    std::vector<double> Delta_fit;
    for (int i = 0; i < n; i++) {
      Delta_fit.push_back(m_D * time_vector[i] + c_D);
    }

    
    // Find Delta_max and Delta_min for Delta_range
    Delta_max = 0;
    Delta_min = 0;
    
    for (int i = 0; i < Delta_fit.size(); i++) {
      if (Delta_fit[i] >= Delta_max) {
        Delta_max = Delta_fit[i];
      }
      if (Delta_fit[i] <= Delta_min) {
        Delta_min = Delta_fit[i];
      }
    }
    
    Delta_range = Delta_max - Delta_min;

    std::cout << "(BEFORE) Delta_max: " << Delta_max << std::endl;
    std::cout << "(BEFORE) Delta_min: " << Delta_min << std::endl;
    std::cout << "(BEFORE) Delta_range: " << Delta_range << std::endl;
  
    SoC1.push_back(next_SoC1);
    SoC2.push_back(next_SoC2);
    Delta.push_back(std::abs(next_SoC1 - next_SoC2));
    std::cout << "Next Delta value: " << std::abs(next_SoC1 - next_SoC2) << std::endl;
    SoC1_predicted = m1 * (time_vector.back()) + c1;
    std::cout << "'SoC1_predicted' variable value: " << SoC1_predicted << std::endl;
    SoC2_predicted = m2 * (time_vector.back()) + c2;
    std::cout << "'SoC2_predicted' variable value: " << SoC2_predicted << std::endl;
    PvA_SoC1 = std::abs(SoC1_predicted - next_SoC1);
    PvA_SoC2 = std::abs(SoC2_predicted - next_SoC2);
    std::cout << "PvA_SoC1: " << PvA_SoC1 << std::endl;
    std::cout << "PvA_SoC2: " << PvA_SoC2 << std::endl;
    
    //////////VALIDITY DETERMINATION////////////////////////-------------------------------------------------
  
    if (Charging == true && Discharging == false){
      if (next_SoC1 < SoC1[SoC1.size()-2]){
        SoC1_INVALID = true;
      }
      else{
        SoC1_INVALID = false;
      }
      if (next_SoC2 < SoC2[SoC2.size()-2]){
        SoC2_INVALID = true;
      }
      else{
        SoC2_INVALID = false;
      }
    }
    if (Charging == false && Discharging == true){
      if (next_SoC1 > SoC1[SoC1.size()-2]){  //-2 because the next_SoCn values have already been appended so not -1
        SoC1_INVALID = true;
      }
      else{
        SoC1_INVALID = false;
      }
      if (next_SoC2 > SoC2[SoC2.size()-2]){
        SoC2_INVALID = true;
      }
      else{
        SoC2_INVALID = false;
      }
    }
    if (PvA_SoC1 >= EPSILON*Delta_range){
      SoC1_INVALID = true;
    }
    else{
        SoC1_INVALID = false;
      }
    if (PvA_SoC2 >= EPSILON*Delta_range){
      SoC2_INVALID = true;
    }
    else{
        SoC2_INVALID = false;
      }


    
    ///////////////////CALIBRATION//////////////////////////////////
  
    

  
    if (SoC1_INVALID == true && SoC2_INVALID == true){  //Trigger Invalid SoC Measurement if all SoC inputs are invalid
      InvalidLoopCounter++;
      next_SoC_calibrated = -50;  //To keep track of invalid SoC_calibrated values, they are -50 with still a time plot - these should be cleaned with the cleanVector function and a separate output .txt file
      if (InvalidLoopCounter == 5){  //Trigger ERRORMSG after 5 Invalid SoC Measurements
        ERRORMSG++;
        std::cout<<"ERROR MESSAGE!!!!! (# " << ERRORMSG << ")" << std::endl;
        if (ERRORMSG == 2){  //Only terminate program after 2 ERRORMSGs (first one is voided to account for decreasing EPSILON hyperparameter)
          outputData (time_vector,SoC1,SoC1_fit,SoC2,SoC2_fit,Delta,Delta_fit,SoC_calibrated,EPSILON_vector,MSE_SoC1,MSE_SoC2,MSE_SoC_calibrated);
          break;
        }
        InvalidLoopCounter = 0;  //Reset Invalid Loop Counter if 1st ERRORMSG is to be voided
      }
    }
    else{
      InvalidLoopCounter = 0;
      if (SoC1_INVALID == false && SoC2_INVALID == true){
        weight1++;
        next_SoC_calibrated = next_SoC1;
      }
      if (SoC1_INVALID == true && SoC2_INVALID == false){
        weight2++;
        next_SoC_calibrated = next_SoC2;
      }
      if (SoC1_INVALID == false && SoC2_INVALID == false){
        weight1++;
        weight2++;
        next_SoC_calibrated = (weight1*next_SoC1 + weight2*next_SoC2)/(weight1+weight2); 
      }
    }



    std::cout << "SoC1_INVALID?: " << SoC1_INVALID << std::endl;
    std::cout << "SoC2_INVALID?: " << SoC2_INVALID << std::endl;
    std::cout << "next_SoC_calibrated: " << next_SoC_calibrated << std::endl;
    
    SoC_calibrated.push_back(next_SoC_calibrated);
    MSE_SoC1.push_back(MSE(SoC1,time_vector));
    MSE_SoC2.push_back(MSE(SoC2,time_vector));
    MSE_SoC_calibrated.push_back(MSE(SoC_calibrated,time_vector));
    

    //End-of-loop increments/decrements
    if (InvalidLoopCounter == 0 && ERRORMSG == 0){  //Only decrement hyperparameter if any of the SoC measurements were valid
      EPSILON = EPSILON - EPSILON_step;
    }
    EPSILON_vector.push_back(EPSILON);
    time_vector.push_back(time_vector.back()+time_step);
    WhileLoopCounter++;
  

  
    std::cout << "EPSILON: " << EPSILON << std::endl;
    std::cout << "InvalidLoopCounter: " << InvalidLoopCounter << std::endl;
    std::cout << "Weight 1: "<< weight1 << std::endl;
    std::cout << "Weight 2: "<< weight2 << std::endl;
    std::cout << "Last time_vector element: " << time_vector.back() << std::endl;



    //Manual iteration cap:
    if (WhileLoopCounter == 10000){
      outputData (time_vector,SoC1,SoC1_fit,SoC2,SoC2_fit,Delta,Delta_fit,SoC_calibrated,EPSILON_vector,MSE_SoC1,MSE_SoC2,MSE_SoC_calibrated);
      break;
    } 
  }
}