
        format none;

        function xdot = f (x, t)
            # Zmena hodnoty x za jednotku casu
            xdot = zeros (9,1);
            xdot(2) = 0; # fixed, kF 
xdot(9) = 0; # fixed, Thylakoid 
xdot(1) = x(1) = x(2) * ( x(5) + x(4) ); # assignment, F(t)  FROM INIT EXPRESSION!!! 
xdot(3) =  + (x(2)*x(5))  + (200000*x(5))  + (x(2)*x(4))  + (200000*x(4)) ; # reaction, P680 / Qa / PhD1 / ac / ChlD1 
xdot(4) =  + (500000000000*x(6))  + (19200000000*x(5)-25000000000*x(4))  - (500000000000*x(4))  - (x(2)*x(4))  - (200000*x(4)) ; # reaction, P680 / Qa / PhD1 / ac / ChlD1* 
xdot(5) =  - (x(2)*x(5))  - (200000*x(5))  - (19200000000*x(5)-25000000000*x(4)) ; # reaction, P680 / Qa / PhD1 / ac* / ChlD1 
xdot(6) =  - (500000000000*x(6))  - (250000000000*x(6))  + (500000000000*x(4))  + (67000000000*x(7)) ; # reaction, P680 / Qa / PhD1- / ac / ChlD1+ 
xdot(7) =  + (250000000000*x(6))  - (4800000000*x(7))  - (67000000000*x(7))  + (2400000000*x(8)) ; # reaction, P680+ / Qa / PhD1- / ac / ChlD1 
xdot(8) =  + (4800000000*x(7))  - (2400000000*x(8)) ; # reaction, P680+ / Qa- / PhD1 / ac / ChlD1 
xdot(1) = 0; # end

        endfunction

        # time scale
        t = logspace(log10(0.0000000000001), log(0.002), 100);

        # default values
        x0 = zeros (9,1);
        x0(2) = 67000000; # fixed parameter, kF
x0(3) = 0; # fixed specie, P680 / Qa / PhD1 / ac / ChlD1
x0(4) = 0; # fixed specie, P680 / Qa / PhD1 / ac / ChlD1*
x0(5) = 1; # fixed specie, P680 / Qa / PhD1 / ac* / ChlD1
x0(6) = 0; # fixed specie, P680 / Qa / PhD1- / ac / ChlD1+
x0(7) = 0; # fixed specie, P680+ / Qa / PhD1- / ac / ChlD1
x0(8) = 0; # fixed specie, P680+ / Qa- / PhD1 / ac / ChlD1
x0(9) = 1; # fixed compartment, Thylakoid
x0(1) = x0(2) * ( x0(5) + x0(4) ); # assignment parameter, F(t)


        #lsode_options ( 'absolute tolerance', max (50*eps, 0.5e-28) );
        lsode_options ( 'relative tolerance', max (50*eps, 0.5e-28) ); 

        y = lsode ("f", x0, t);
        
        for i = 1:length(t)
            y(i, 1) = y(i, 2) * ( y(i, 5) + y(i, 4) ); # assignment, F(t)  FROM INIT EXPRESSION!!! 

        endfor
        
        [t',y]
        