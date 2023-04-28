
clear
% param
Par=[1150; 0.0012; 0.0012; 0.0023; 0.0023; 3300 * 5.5;  6.9 * 5.68; 6.2 * 5.68; 0.017 * 5.68; 0.011 * 5.65; 0.4005; 0.0012; 0.0017;  0.0012;  0.0023; 0.0023;  0.0012;   0.0012;     1.9E-7;      4.7E-7; 45 * 5.68;   36.7 ;    0.0012;     0.0023;    0.0012;    0.0023; 6.16; 0.000356 ;     0.7;   0.000009088;  0.006;  0.6; 1; 1; 0.5;  0.5; 0.5; 0.5;0.5; -1; 0.5; 0.5; 0.5; 0.5; 0.5];  

% initial conditions 
x0 = [100  100  100  100  100 100 100 100 100 ...
      29.08138  88.89034   70.88889 73.48281 76.90346 ...
     10993.27 5608.411];

tic
tspan = 1:10000; 
[t,x] = ode23s(@CholDifEqlf,tspan,x0',[],Par); 
toc %26s

figure(1)
plot(t,x(:,15))
xlabel('time')
ylabel('Hepatic Cholesterol')



function dxdt = CholDifEqlf(t,x,p)
%param
r9 = p(1); 
r1 = p(2);
r7 = p(3);
r11 = p(4);  
r17 = p(5); 
r19 = p(6); 
r20 = p(7) ; 	
r21= p(8); 	
r23 = p(9); 	  	 
r24 = p(10); 
r25= p(11); 
r3=p(12) ;    
r13= p(13);   					
r5 = p(14); 				
r15 = p(15); 
r10= p(16);  
r2= p(17); 		  
r8= p(18); 	 
r12= p(19); 	
r18 = p(20); 				
r22 = p(21); 							 										
r28 = p(22); 							
r4 = p(23); 	
r14 = p(24);
r6 = p(25);
r16 = p(26);  	
r26 = p(27);   
r27 = p(28);         
r30 = p(29);   			
r33 = p(30);    
r31 = p(31);  
r32 = p(32); 
DE = p(33); 
DEP = p(34);
a5 = p(35); 
a6 = p(36);
a7 = p(37); 
a8 = p(38); 
a10 = p(39); 
a22 = p(40); 
a24 = p(41); 
a25 = p(42); 
a29 = p(43); 
a32 = p(44); 
aPL = p(45); 


dxdt = zeros (16,1);
%HMGCR mRNA
dxdt(1) = r1 * x(5)  - r2 * x(1);
%FDPS mRNA
dxdt(2) = r3 * x(5)  - r4 * x(2); 
%Lanosterol synthase mRNA
dxdt(3) = r5* x(5)  - r6 * x(3);
%Squalene Synthase mRNA 
dxdt(4) = r7 * x(5)  - r8 * x(4) ; 
%Proteins
%SREBP
dxdt(5)  = r9 * x(15) ^ a22 - r10 * x(5) ;
%HMG CoA Reductase
dxdt(6)  = r11 * x(1)  - r12 * x(6)  * x(11) ^ a5 * x(12) ^ a6 * x(14) ^ a7 * x(15) ^ a8 ;
% Farnesyl diphosphate synthase 
dxdt(7)  = r13 * x(2)  - r14 * x(7) ;
%Lanosterol Synthase 
dxdt(8)  = r15 * x(3)  - r16 * x(8) ;
%Squalene Synthase
dxdt(9)  = r17 * x(5)  - r18 * x(9) * x(15);
%Metabolites
%HMGCoA
dxdt(10)  = (r19 - r20 *  x(10)  * x(6)^ a10 *  x(7) ^ a25) * (1000/23);
%Geranyl-PP
dxdt(11)  =  (r20 *  x(10)   * x(6) ^ a10 * x(7) ^ a25 - r21 * x(11)  * x(7) ^ a29 - r33 * x(11)) * (1000/23);
%Farnesyl- PP
dxdt(12)  =  (r21 * x(11) * x(7) ^ a29  - r22 * x(12) - r23 * x(12) * x(9) ^ a24) * (1000/23);
%squalene
dxdt(13)  = (r23 * x(12) * x(9) ^ a24 - r24 * x(13) * x(8) ^ a32) * (1000/23);
%lanesterol
dxdt(14)  = (r24  * x(13) * x(8) ^ a32  - r25 *  x(14) ) * (1000/23) ;
%cholesterol
dxdt(15)  =   (r25 * x(14)  + r26 - r27 * x(15) - r28 * x(15) ^ aPL  +   r30 * x(16) ) * (1000/23);
dxdt(16) = r28 * x(15) ^ aPL - r30 * x(16) - r31 * x(16) + r32;
end