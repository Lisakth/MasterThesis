
three_compartment_5k_fengconv_stan <- "

real Integral_Feng_3TC5K(real time, real A, real B, real C, 
                       real alpha, real beta, real gamma, 
                       real Th1, real Th2, real Ph1, real Ph2, real Ph3) {
     
    real integ;
    
    integ = (B*(1 - exp(-(alpha*time)))*Ph3)/alpha^2 + (C*(1 - exp(-(alpha*time)))*Ph3)/alpha^2 + (A*(-1 + exp(-(alpha*time)))*Ph3)/alpha^3 + 
  (B*(-1 + exp(-(beta*time)))*Ph3)/beta^2 + (C*(-1 + exp(-(gamma*time)))*Ph3)/gamma^2 + (A*Ph3*time)/alpha^2 - (B*Ph3*time)/alpha + (B*Ph3*time)/beta - 
  (C*Ph3*time)/alpha + (C*Ph3*time)/gamma - (A*Ph3*(1 - (1 + alpha*time)/exp(alpha*time)))/alpha^3 - (A*(1 - exp(-(alpha*time)))*Ph1)/(alpha*(alpha-Th1)^2) - 
  (A*Ph1*(1 - (1 + alpha*time)/exp(alpha*time)))/(alpha*(alpha-Th1)^2) + (B*(1 - exp(-(alpha*time)))*Ph1)/(alpha*(alpha - Th1)) + 
  (C*(1 - exp(-(alpha*time)))*Ph1)/(alpha*(alpha - Th1)) + (A*(1 - exp(-(time*Th1)))*Ph1)/((alpha-Th1)^2*Th1) - 
  (B*(1 - exp(-(time*Th1)))*Ph1)/((alpha - Th1)*Th1) - (C*(1 - exp(-(time*Th1)))*Ph1)/((alpha - Th1)*Th1) + (B*(1 - exp(-(time*Th1)))*Ph1)/((beta - Th1)*Th1) + 
  (C*(1 - exp(-(time*Th1)))*Ph1)/((gamma - Th1)*Th1) + (A*Ph1*(1 - (1 + alpha*time)/exp(alpha*time))*Th1)/(alpha^2*(alpha-Th1)^2) + 
  (B*(1 - exp(-(beta*time)))*Ph1)/(beta*(-beta + Th1)) + (C*(1 - exp(-(gamma*time)))*Ph1)/(gamma*(-gamma + Th1)) - 
  (A*(1 - exp(-(alpha*time)))*Ph2)/(alpha*(alpha-Th2)^2) - (A*Ph2*(1 - (1 + alpha*time)/exp(alpha*time)))/(alpha*(alpha-Th2)^2) + 
  (B*(1 - exp(-(alpha*time)))*Ph2)/(alpha*(alpha - Th2)) + (C*(1 - exp(-(alpha*time)))*Ph2)/(alpha*(alpha - Th2)) + 
  (A*(1 - exp(-(time*Th2)))*Ph2)/((alpha-Th2)^2*Th2) - (B*(1 - exp(-(time*Th2)))*Ph2)/((alpha - Th2)*Th2) - (C*(1 - exp(-(time*Th2)))*Ph2)/((alpha - Th2)*Th2) + 
  (B*(1 - exp(-(time*Th2)))*Ph2)/((beta - Th2)*Th2) + (C*(1 - exp(-(time*Th2)))*Ph2)/((gamma - Th2)*Th2) + 
  (A*Ph2*(1 - (1 + alpha*time)/exp(alpha*time))*Th2)/(alpha^2*(alpha-Th2)^2) + (B*(1 - exp(-(beta*time)))*Ph2)/(beta*(-beta + Th2)) + 
  (C*(1 - exp(-(gamma*time)))*Ph2)/(gamma*(-gamma + Th2));
   
   return(integ);
                       
}


real threetc5kND_log_stan_BPp(real logK1, real logVnd, real logBPp, 
                    real logk4, real logk5, real logvB, real time,
                    real t0, 
                    real A, real B, real C, 
                    real alpha, real beta, real gamma,
                    real ti, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real k5;
    real KI;
    real Vnd;
    real BPnd;
    real BPp;
    real vB;
    
    real Th1;
    real Th2;
    real Ph1;
    real Ph2;
    real Ph3;
    
    real tcorr;
    
    real pred;
    

    
    K1 = exp(logK1);
    Vnd = exp(logVnd);
    BPp = exp(logBPp);
    k4 = exp(logk4);
    vB = exp(logvB);
    k5 = exp(logk5);
    
    k2 = K1 / Vnd - k5;
    k3 = (k4 * (k2 + k5) * BPp)/K1;
    KI = K1 * k5 / (k2 + k5);

    Th1 = 0.5 * (k2 + k3 + k4 + k5 - sqrt((k2 + k3 + k4 + k5)^2 - 4 * k4 * (k2+k5)));
    Th2 = 0.5 * (k2 + k3 + k4 + k5 + sqrt((k2 + k3 + k4 + k5)^2 - 4 * k4 * (k2+k5)));

    Ph1 = (K1 * (k3 + k4 + k5 - Th1) - KI * Th2)/(Th2 - Th1);
    Ph2 = (K1 * (Th2 - k3 - k4 - k5) + KI * Th1)/(Th2 - Th1);
    Ph3 = KI;
    
    tcorr = time - t0;
  
  pred =  // 2TC contribution
          (1-vB) * (
            (tcorr > 0) * (             // Before t0 equal to 0
              
              // Before t0+ti
              (tcorr < ti)*(1/ti)*(     
                Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3)
              ) +
                
              // After t0+ti 
              (tcorr >= ti)*(1/ti)* (   
                Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3) -
                  Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3)
              )
            )
          ) +
          // Blood contribution
          vB*bloodval;
       
  return(pred);
}

// S

real threetc5kS_log_stan_BPp(real logK1, real logVnd, real logBPp, 
                    real logk4, real logk5, real logvB, real time,
                    real t0, 
                    real A, real B, real C, 
                    real alpha, real beta, real gamma,
                    real ti, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real k5;
    real KI;
    real Vnd;
    real BPnd;
    real BPp;
    real vB;
    
    real Th1;
    real Th2;
    real Ph1;
    real Ph2;
    real Ph3;
    
    real tcorr;
    
    real pred;
    

    
    K1 = exp(logK1);
    Vnd = exp(logVnd);
    BPp = exp(logBPp);
    k4 = exp(logk4);
    vB = exp(logvB);
    k5 = exp(logk5);
    
    k2 = (K1 - k5 * BPp)/Vnd;
    k3 = BPp * (k4 + k5)/Vnd;
    KI = K1 * k3 * k5 / (k2 * k4 + k2 * k5 + k3 * k5);
  
    Th1 = 0.5 * (k2 + k3 + k4 + k5 + sqrt((k2 + k3 + k4 + k5)^2 - 4 * (k2 * k4 + k2 * k5 + k3 * k5)));
    Th2 = 0.5 * (k2 + k3 + k4 + k5 - sqrt((k2 + k3 + k4 + k5)^2 - 4 * (k2 * k4 + k2 * k5 + k3 * k5)));

    Ph1 = (K1 * (k3 + k4 + k5 - Th1) - KI * Th2)/(Th2 - Th1);
    Ph2 = (K1 * (Th2 - k3 - k4 - k5) + KI * Th1)/(Th2 - Th1);
    Ph3 = KI;
    
    tcorr = time - t0;
  
  pred =  // 2TC contribution
          (1-vB) * (
            (tcorr > 0) * (             // Before t0 equal to 0
              
              // Before t0+ti
              (tcorr < ti)*(1/ti)*(     
                Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3)
              ) +
                
              // After t0+ti 
              (tcorr >= ti)*(1/ti)* (   
                Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3) -
                  Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3)
              )
            )
          ) +
          // Blood contribution
          vB*bloodval;
       
  return(pred);
}


// P

real threetc5kP_log_stan_BPp(real logK1, real logVnd, real logBPp, 
                    real logk4, real logk5, real logvB, real time,
                    real t0, 
                    real A, real B, real C, 
                    real alpha, real beta, real gamma,
                    real ti, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real k5;
    real KI;
    real Vnd;
    real BPnd;
    real BPp;
    real vB;
    
    real Th1;
    real Th2;
    real Ph1;
    real Ph2;
    real Ph3;
    
    real tcorr;
    
    real pred;
    

    
    K1 = exp(logK1);
    Vnd = exp(logVnd);
    BPp = exp(logBPp);
    k4 = exp(logk4);
    vB = exp(logvB);
    k5 = exp(logk5);
    
    k2 = K1 / Vnd;
    k3 = (k4 * BPp)/Vnd;
    KI = k5;
  
    Th1 = 0.5 * (k2 + k3 + k4 - sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4));
    Th2 = 0.5 * (k2 + k3 + k4 + sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4));

    Ph1 = (K1 * (k3 + k4 - Th1))/(Th2 - Th1);
    Ph2 = (K1 * (Th2 - k3 - k4))/(Th2 - Th1);
    Ph3 = KI;
    
    tcorr = time - t0;
  
  pred =  // 2TC contribution
          (1-vB) * (
            (tcorr > 0) * (             // Before t0 equal to 0
              
              // Before t0+ti
              (tcorr < ti)*(1/ti)*(     
                Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3)
              ) +
                
              // After t0+ti 
              (tcorr >= ti)*(1/ti)* (   
                Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3) -
                  Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma, 
                                 Th1, Th2, Ph1, Ph2, Ph3)
              )
            )
          ) +
          // Blood contribution
          vB*bloodval;
       
  return(pred);
}
                    
//
//
//
//
//
//
//
//
//
//
//
//
//
// mixture models
//
//
//
//
//
//
//
//
//
//
//
//
//
//
real threetc5kND_log_stan_mix_BPp(real logK1, real logVnd, real logBPp, real logk4,
                    real logk5, real logvB,
                    real logBPpdiff, real logk4diff, real logk5diff,
                    real prophigh, real proplow,
                    real time,
                    real t0,
                    real A, real B, real C,
                    real alpha, real beta, real gamma,
                    real ti, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real k5;
    real KI;
    real Vnd;
    real BPp;
    real vB;

    real Th1;
    real Th2;
    real Ph1;
    real Ph2;
    real Ph3;
    
    real Th1low;
    real Th2low;
    real Ph1low;
    real Ph2low;
    real Ph3low;

    real tcorr;

    real predhigh;
    real predlow;
    real pred;
    
    real BPplow;
    real k4low;
    real k3low;
    real k5low;
    real KIlow;
    real k5mean;
    

    K1 = exp(logK1);
    Vnd = exp(logVnd);
    vB = exp(logvB);

    // High
    BPp = exp(logBPp);
    k4 = exp(logk4);
    k5 = exp(logk5);
    
    // Low
    BPplow = exp(logBPp + logBPpdiff);
    k4low = exp(logk4 + logk4diff);
    k5low = exp(logk5 + logk5diff);
    
    // Assumed to be indepenpendent of genotypes
    // note for myself: same no mather if mean(k5) or mean(k2)
    k5mean = (k5 +k5low)/2;
    k2 = K1 / Vnd - k5mean;

    // High
    k3 = (k4 * (k2 + k5) * BPp)/K1;
    KI = K1 * k5 / (k2 + k5);

    Th1 = 0.5 * (k2 + k3 + k4 + k5 - sqrt((k2 + k3 + k4 + k5)^2 - 4 * k4 * (k2+k5)));
    Th2 = 0.5 * (k2 + k3 + k4 + k5 + sqrt((k2 + k3 + k4 + k5)^2 - 4 * k4 * (k2+k5)));

    Ph1 = (K1 * (k3 + k4 + k5 - Th1) - KI * Th2)/(Th2 - Th1);
    Ph2 = (K1 * (Th2 - k3 - k4 - k5) + KI * Th1)/(Th2 - Th1);
    Ph3 = KI;
    
    
    // Low
    k3low = (k4low * (k2 + k5low) * BPplow)/K1;
    KIlow = K1 * k5low / (k2 + k5low);
    
    Th1low = 0.5 * (k2 + k3low + k4low + k5low - sqrt((k2 + k3low + k4low + k5low)^2 - 4 * k4low * (k2 + k5low)));
    Th2low = 0.5 * (k2 + k3low + k4low + k5low + sqrt((k2 + k3low + k4low + k5low)^2 - 4 * k4low * (k2 + k5low)));

    Ph1low = (K1 * (k3low + k4low + k5low - Th1low) - KIlow * Th2low)/(Th2low - Th1low);
    Ph2low = (K1 * (Th2low - k3low - k4low - k5low) + KIlow * Th1low)/(Th2low - Th1low);
    Ph3low = KIlow;
    

    tcorr = time - t0;

  predhigh =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC off CND
      
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2,Ph3)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3)
        )
      ) 
      
  ) +
  // Blood contribution
  vB *  bloodval;
  
  predlow =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        )
      )
  ) +
  // Blood contribution
  vB *  bloodval;

  pred = ( prophigh * predhigh ) + ( proplow * predlow );

  return(pred);
}
                    
real threetc5kS_log_stan_mix_BPp(real logK1, real logVnd, real logBPp, real logk4,
                    real logk5, real logvB,
                    real logBPpdiff, real logk4diff, real logk5diff,
                    real prophigh, real proplow,
                    real time,
                    real t0,
                    real A, real B, real C,
                    real alpha, real beta, real gamma,
                    real ti, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real k5;
    real KI;
    real Vnd;
    real BPp;
    real vB;

    real Th1;
    real Th2;
    real Ph1;
    real Ph2;
    real Ph3;
    
    real Th1low;
    real Th2low;
    real Ph1low;
    real Ph2low;
    real Ph3low;

    real tcorr;

    real predhigh;
    real predlow;
    real pred;
    
    real BPplow;
    real k4low;
    real k3low;
    real k5low;
    real KIlow;
    
    real k5mean;
    real BPpmean;
    

    K1 = exp(logK1);
    Vnd = exp(logVnd);
    vB = exp(logvB);

    // High
    BPp = exp(logBPp);
    k4 = exp(logk4);
    k5 = exp(logk5);
    
    // Low
    BPplow = exp(logBPp + logBPpdiff);
    k4low = exp(logk4 + logk4diff);
    k5low = exp(logk5 + logk5diff);
    
    // Assumed to be indepenpendent of genotypes
    // note for myself: same no mather if mean(k5) or mean(k2)
    k5mean = (k5 + k5low)/2;
    BPpmean = (BPp + BPplow)/2;

    k2 = (K1 - k5mean * BPpmean)/Vnd;

    
    // High

    k3 = BPp * (k4 + k5)/Vnd;
    KI = K1 * k3 * k5 / (k2 * k4 + k2 * k5 + k3 * k5);
  
    Th1 = 0.5 * (k2 + k3 + k4 + k5 + sqrt((k2 + k3 + k4 + k5)^2 - 4 * (k2 * k4 + k2 * k5 + k3 * k5)));
    Th2 = 0.5 * (k2 + k3 + k4 + k5 - sqrt((k2 + k3 + k4 + k5)^2 - 4 * (k2 * k4 + k2 * k5 + k3 * k5)));

    Ph1 = (K1 * (k3 + k4 + k5 - Th1) - KI * Th2)/(Th2 - Th1);
    Ph2 = (K1 * (Th2 - k3 - k4 - k5) + KI * Th1)/(Th2 - Th1);
    Ph3 = KI;
    
    
    // Low
    
    k3low = BPplow * (k4low + k5low)/Vnd;
    KIlow = K1 * k3low * k5low / (k2 * k4low + k2 * k5low + k3low * k5low);
    
    Th1low = 0.5 * (k2 + k3low + k4low + k5low + sqrt((k2 + k3low + k4low + k5low)^2 - 4 * (k2 * k4low + k2 * k5low + k3low * k5low)));
    Th2low = 0.5 * (k2 + k3low + k4low + k5low - sqrt((k2 + k3low + k4low + k5low)^2 - 4 * (k2 * k4low + k2 * k5low + k3low * k5low)));

    Ph1low = (K1 * (k3low + k4low + k5low - Th1low) - KIlow * Th2low)/(Th2low - Th1low);
    Ph2low = (K1 * (Th2low - k3low - k4low - k5low) + KIlow * Th1low)/(Th2low - Th1low);
    Ph3low = KIlow;

    

    tcorr = time - t0;

  predhigh =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC off CND
      
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2,Ph3)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3)
        )
      ) 
      
  ) +
  // Blood contribution
  vB *  bloodval;
  
  predlow =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        )
      )
  ) +
  // Blood contribution
  vB *  bloodval;

  pred = ( prophigh * predhigh ) + ( proplow * predlow );

  return(pred);
}

real threetc5kP_log_stan_mix_BPp(real logK1, real logVnd, real logBPp, real logk4,
                    real logk5, real logvB,
                    real logBPpdiff, real logk4diff, real logk5diff,
                    real prophigh, real proplow,
                    real time,
                    real t0,
                    real A, real B, real C,
                    real alpha, real beta, real gamma,
                    real ti, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real k5;
    real KI;
    real Vnd;
    real BPp;
    real vB;

    real Th1;
    real Th2;
    real Ph1;
    real Ph2;
    real Ph3;
    
    real Th1low;
    real Th2low;
    real Ph1low;
    real Ph2low;
    real Ph3low;

    real tcorr;

    real predhigh;
    real predlow;
    real pred;
    
    real BPplow;
    real k4low;
    real k3low;
    real k5low;
    real KIlow;
    

    K1 = exp(logK1);
    Vnd = exp(logVnd);
    vB = exp(logvB);

    // High
    BPp = exp(logBPp);
    k4 = exp(logk4);
    k5 = exp(logk5);


    k2 = K1 / Vnd;
    k3 = (k4 * BPp)/Vnd;
    KI = k5;
  
    Th1 = 0.5 * (k2 + k3 + k4 - sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4));
    Th2 = 0.5 * (k2 + k3 + k4 + sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4));

    Ph1 = (K1 * (k3 + k4 - Th1))/(Th2 - Th1);
    Ph2 = (K1 * (Th2 - k3 - k4))/(Th2 - Th1);
    Ph3 = KI;
    
    
    // Low
    BPplow = exp(logBPp + logBPpdiff);
    k4low = exp(logk4 + logk4diff);
    k5low = exp(logk5 + logk5diff);
    
    k3low = (k4low * BPplow)/Vnd;
    KIlow = k5low;
    
    Th1low = 0.5 * (k2 + k3low + k4low - sqrt((k2 + k3low + k4low)^2 - 4 * k4low * k2));
    Th2low = 0.5 * (k2 + k3low + k4low + sqrt((k2 + k3low + k4low)^2 - 4 * k4low * k2));

    Ph1low = (K1 * (k3low + k4low - Th1low))/(Th2low - Th1low);
    Ph2low = (K1 * (Th2low - k3low - k4low))/(Th2low - Th1low);
    Ph3low = KIlow;

    

    tcorr = time - t0;

  predhigh =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC off CND
      
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2,Ph3)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3)
        )
      ) 
      
  ) +
  // Blood contribution
  vB *  bloodval;
  
  predlow =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        )
      )
  ) +
  // Blood contribution
  vB *  bloodval;

  pred = ( prophigh * predhigh ) + ( proplow * predlow );

  return(pred);
}

real threetc5kND_log_stan_mix_BPp_k2(real logK1, real logVnd, real logBPp, real logk4,
                    real logk5, real logvB,
                    real logBPpdiff, real logk4diff, real logk5diff, real logVnddiff,
                    real prophigh, real proplow,
                    real time,
                    real t0,
                    real A, real B, real C,
                    real alpha, real beta, real gamma,
                    real ti, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real k5;
    real KI;
    real Vnd;
    real BPp;
    real vB;

    real Th1;
    real Th2;
    real Ph1;
    real Ph2;
    real Ph3;
    
    real Th1low;
    real Th2low;
    real Ph1low;
    real Ph2low;
    real Ph3low;

    real tcorr;

    real predhigh;
    real predlow;
    real pred;
    
    real BPplow;
    real k4low;
    real k3low;
    real k5low;
    real KIlow;
    real k5mean;
    real Vndlow;
    real k2low;
    

    K1 = exp(logK1);
    vB = exp(logvB);

    // High
    BPp = exp(logBPp);
    k4 = exp(logk4);
    k5 = exp(logk5);
    Vnd = exp(logVnd);
    
    // Low
    BPplow = exp(logBPp + logBPpdiff);
    k4low = exp(logk4 + logk4diff);
    k5low = exp(logk5 + logk5diff);
    Vndlow = exp(logVnd + logVnddiff);
    
    // Assumed to be indepenpendent of genotypes
    // note for myself: same no mather if mean(k5) or mean(k2)
    // k5mean = (k5 +k5low)/2;
    // k2 = K1 / Vnd - k5mean;

    // High
    k2 = K1 / Vnd - k5;
    k3 = (k4 * (k2 + k5) * BPp)/K1;
    KI = K1 * k5 / (k2 + k5);

    Th1 = 0.5 * (k2 + k3 + k4 + k5 - sqrt((k2 + k3 + k4 + k5)^2 - 4 * k4 * (k2+k5)));
    Th2 = 0.5 * (k2 + k3 + k4 + k5 + sqrt((k2 + k3 + k4 + k5)^2 - 4 * k4 * (k2+k5)));

    Ph1 = (K1 * (k3 + k4 + k5 - Th1) - KI * Th2)/(Th2 - Th1);
    Ph2 = (K1 * (Th2 - k3 - k4 - k5) + KI * Th1)/(Th2 - Th1);
    Ph3 = KI;
    
    
    // Low
    k2low = K1 / Vndlow - k5low;
    k3low = (k4low * (k2low + k5low) * BPplow)/K1;
    KIlow = K1 * k5low / (k2low + k5low);
    
    Th1low = 0.5 * (k2low + k3low + k4low + k5low - sqrt((k2low + k3low + k4low + k5low)^2 - 4 * k4low * (k2low + k5low)));
    Th2low = 0.5 * (k2low + k3low + k4low + k5low + sqrt((k2low + k3low + k4low + k5low)^2 - 4 * k4low * (k2low + k5low)));

    Ph1low = (K1 * (k3low + k4low + k5low - Th1low) - KIlow * Th2low)/(Th2low - Th1low);
    Ph2low = (K1 * (Th2low - k3low - k4low - k5low) + KIlow * Th1low)/(Th2low - Th1low);
    Ph3low = KIlow;
    

    tcorr = time - t0;

  predhigh =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC off CND
      
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2,Ph3)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3)
        )
      ) 
      
  ) +
  // Blood contribution
  vB *  bloodval;
  
  predlow =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        )
      )
  ) +
  // Blood contribution
  vB *  bloodval;

  pred = ( prophigh * predhigh ) + ( proplow * predlow );

  return(pred);
}
                    
real threetc5kS_log_stan_mix_BPp_k2(real logK1, real logVnd, real logBPp, real logk4,
                    real logk5, real logvB,
                    real logBPpdiff, real logk4diff, real logk5diff, real logVnddiff,
                    real prophigh, real proplow,
                    real time,
                    real t0,
                    real A, real B, real C,
                    real alpha, real beta, real gamma,
                    real ti, real bloodval) {

    real K1;
    real k2;
    real k3;
    real k4;
    real k5;
    real KI;
    real Vnd;
    real BPp;
    real vB;

    real Th1;
    real Th2;
    real Ph1;
    real Ph2;
    real Ph3;
    
    real Th1low;
    real Th2low;
    real Ph1low;
    real Ph2low;
    real Ph3low;

    real tcorr;

    real predhigh;
    real predlow;
    real pred;
    
    real BPplow;
    real k4low;
    real k3low;
    real k5low;
    real KIlow;
    real Vndlow;
    real k2low;
    
    // real k5mean;
    // real BPpmean;
    

    K1 = exp(logK1);
    vB = exp(logvB);

    // High
    BPp = exp(logBPp);
    k4 = exp(logk4);
    k5 = exp(logk5);
    Vnd = exp(logVnd);
    
    // Low
    BPplow = exp(logBPp + logBPpdiff);
    k4low = exp(logk4 + logk4diff);
    k5low = exp(logk5 + logk5diff);
    Vndlow = exp(logVnd + logVnddiff);
    
    // Assumed to be indepenpendent of genotypes
    // note for myself: same no mather if mean(k5) or mean(k2)
    // k5mean = (k5 + k5low)/2;
    // BPpmean = (BPp + BPplow)/2;

    
    // High

    k2 = (K1 - k5 * BPp)/Vnd;
    k3 = BPp * (k4 + k5)/Vnd;
    KI = K1 * k3 * k5 / (k2 * k4 + k2 * k5 + k3 * k5);
  
    Th1 = 0.5 * (k2 + k3 + k4 + k5 + sqrt((k2 + k3 + k4 + k5)^2 - 4 * (k2 * k4 + k2 * k5 + k3 * k5)));
    Th2 = 0.5 * (k2 + k3 + k4 + k5 - sqrt((k2 + k3 + k4 + k5)^2 - 4 * (k2 * k4 + k2 * k5 + k3 * k5)));

    Ph1 = (K1 * (k3 + k4 + k5 - Th1) - KI * Th2)/(Th2 - Th1);
    Ph2 = (K1 * (Th2 - k3 - k4 - k5) + KI * Th1)/(Th2 - Th1);
    Ph3 = KI;
    
    
    // Low
    
    k2low = (K1 - k5low * BPplow)/Vndlow;
    k3low = BPplow * (k4low + k5low)/Vndlow;
    KIlow = K1 * k3low * k5low / (k2low * k4low + k2low * k5low + k3low * k5low);
    
    Th1low = 0.5 * (k2low + k3low + k4low + k5low + sqrt((k2low + k3low + k4low + k5low)^2 - 4 * (k2low * k4low + k2low * k5low + k3low * k5low)));
    Th2low = 0.5 * (k2low + k3low + k4low + k5low - sqrt((k2low + k3low + k4low + k5low)^2 - 4 * (k2low * k4low + k2low * k5low + k3low * k5low)));

    Ph1low = (K1 * (k3low + k4low + k5low - Th1low) - KIlow * Th2low)/(Th2low - Th1low);
    Ph2low = (K1 * (Th2low - k3low - k4low - k5low) + KIlow * Th1low)/(Th2low - Th1low);
    Ph3low = KIlow;

    

    tcorr = time - t0;

  predhigh =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC off CND
      
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2,Ph3)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1, Th2, Ph1, Ph2, Ph3)
        )
      ) 
      
  ) +
  // Blood contribution
  vB *  bloodval;
  
  predlow =  
  (1-vB) * (
    (tcorr > 0) * (             // Before t0 equal to 0

      // 2TC
        // Before t0+ti
        (tcorr < ti)*(1/ti)*(
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        ) +

        // After t0+ti
        (tcorr >= ti)*(1/ti)* (
          Integral_Feng_3TC5K(tcorr, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low) -
            Integral_Feng_3TC5K(tcorr - ti, A, B, C, alpha, beta, gamma,
                           Th1low, Th2low, Ph1low, Ph2low, Ph3low)
        )
      )
  ) +
  // Blood contribution
  vB *  bloodval;

  pred = ( prophigh * predhigh ) + ( proplow * predlow );

  return(pred);
}

"
