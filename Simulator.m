function [x_out] = Simulator(simPara)
% a stochastic diffiential equation integrator
% author: Yuguang Yang yyang60@jhu.edu

    x = simPara.x0;
    x_out = [x];
    for step = 1: simPara.nstep
        simPara.step = step;
        
        diffusivity_coeff = simPara.diffusivity(x,simPara);
        
        u = diffusivity_coeff*simPara.mobility * simPara.force(x,simPara);
        
        random_disp = chol(2.0*simPara.D*diffusivity_coeff*simPara.dt,'lower')*randn(3,1);
        x = x + u*simPara.dt + random_disp;
        if (mod(step,simPara.saveInterval) == 0 )
            x_out = [x_out x];
        end
    end

