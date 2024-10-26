function logKv = logBesselK(nu, z)
%     % Initial value nu0 is chosen such that Kν0(z) can be computed with high accuracy
%     nu0 = floor(nu); % This might be adjusted based on the stability region.
%     if nu0 == nu
%         logKv = log(besselk(nu, z));
%         return;
%     end
nu0=10;
nu_length = length(nu);
if nu_length == 1
    temp = log(besselk(nu,z));
    if ~isinf(temp)
        logKv = temp;
    else
        while besselk(nu0,z) == 0
            nu0 = nu0 + 10;
        end
        while isinf(besselk(nu0,z))
            nu0 = nu0 - 10;
        end
        % Compute the initial values using MATLAB's built-in function
        K_nu0 = besselk(nu0, z);
        %K_nu0p1 = besselk(nu0+1, z);

        % Use log to avoid underflow/overflow
        logK_nu0 = log(K_nu0);
        %logK_nu0p1 = log(K_nu0p1);
        r_prev = besselk(nu0+1, z) / besselk(nu0, z);
        % Summing the logarithm of ratios using forward recursion
        logKv = logK_nu0;
        for k = 0:(nu-nu0-1)
            r = 1/r_prev+2*(nu0+k)/z;
            log_r = log(r);
            logKv = logKv + log_r;

            % Update values for next iteration
            r_prev = r;
        end
    end
else
    logKv = zeros(nu_length,1);
    for ii = 1:nu_length
        nu0=10;
        nui = nu(ii);zi = z(ii);
        temp = log(besselk(nui,zi));
        if ~isinf(temp)
            logKv(ii) = temp;
        else
            while besselk(nu0,zi) == 0
                nu0 = nu0 + 10;
            end
            while isinf(besselk(nu0,zi))
                nu0 = nu0 - 10;
            end
            % Compute the initial values using MATLAB's built-in function
            K_nu0 = besselk(nu0, zi);
            %K_nu0p1 = besselk(nu0+1, zi);

            % Use log to avoid underflow/overflow
            logK_nu0 = log(K_nu0);
            %logK_nu0p1 = log(K_nu0p1);
            r_prev = besselk(nu0+1, zi) / besselk(nu0, zi);
            % Summing the logarithm of ratios using forward recursion
            logKv(ii) = logK_nu0;
            for k = 0:(nui-nu0-1)
                r = 1/r_prev+2*nui/zi;
                log_r = log(r);
                logKv(ii) = logKv(ii) + log_r;

                % Update values for next iteration
                r_prev = r;
            end
        end
    end
end

end
