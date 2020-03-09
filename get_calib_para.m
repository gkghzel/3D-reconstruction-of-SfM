function [Calib_P1, Calib_P2, FrmNum] = get_calib_para(setName, rightMain) %#ok<INUSD>

    rightMain = (nargin==2);
    if(strcmp(setName,'set1')||strcmp(setName,'set2'))||strcmp(setName, 'set3')
        % intrisic paras
        fc1 = [530.9002; 581.00362]; % Focal_Length
        cc1 = [136.63037; 161.32884] + 1; % Principal_Point
        kc1 = [-0.2865; 0.29524; -0.00212; 0.00152; 0]; % Distortion_coefficient
        alpha_c1 = 0;
        fc2 = [524.84413; 577.11024];
        cc2 = [216.17358; 149.76379] + 1;
        kc2 = [-0.25745; 0.62307; 0.0366; -0.01082; 0];
        alpha_c2 = 0;
        % extrinsic paras 
        om2 = [-0.009954307627539; -0.042566273591172; 0.011454074192639]; % rotation vector right camera to left
        R2 = rodrigues(om2); % rotation matrix
        T2 = [-5.49238; 0.04267; -0.39886]; % transfer para right camera to left
    elseif (strcmp(setName,'set4'))||(strcmp(setName,'set5')) % set4=f5; set5=f7
        % intrisic paras
        fc1 = [391.656525; 426.835144]; % Focal_Length
        cc1 = [165.964371; 154.498138] + 1; % Principal_Point
        kc1 = [-0.196312; 0.129540; 0.004356; 0.006236; 0]; % Distortion_coefficient
        alpha_c1 = 0;
        fc2 = [390.376862; 426.228882]; % Focal_Length
        cc2 = [190.896454; 145.071411] + 1; % Principal_Point
        kc2 = [-0.205824; 0.186125; 0.015374; 0.003660; 0]; % Distortion_coefficient
        alpha_c2 = 0;
        % extrinsic paras 
        R2 = [0.999999 -0.001045 -0.000000;
              0.001045 0.999999  -0.000000;
              0.000000 0.000000  1.000000]; % rotation matrix
        T2 = [-5.520739; -0.031516; -0.051285]; % transfer para right camera to left
    else
        disp(['error: can not find dataset: ', setName]);
        Calib_P1 = 0;
        Calib_P2 = 0;
        return;
    end

    A1 = [fc1(1)    0           cc1(1);
          0         fc1(2)      cc1(2); 
          0         0           1];
    
    A2 = [fc2(1)    0       cc2(1); 
          0         fc2(2), cc2(2); 
          0         0       1];
    
    
    if ~rightMain
        C1 = [A1 zeros(3,1)];
        C2 = A2*[R2, T2]; % 3x4 calibration matrix for right camera A*[R,T], where om is rotation vector and rodrigues() transform it into rotation matrix 
        Calib_P1 = struct('focus', fc1, 'center', cc1, 'distort', kc1, 'alpha_c', alpha_c1, 'matrix', C1);
        Calib_P2 = struct('focus', fc2, 'center', cc2, 'distort', kc2, 'alpha_c', alpha_c2, 'matrix', C2);
    else
        C2 = [A2 zeros(3,1)];
        C1 = A1*[inv(R2), -T2]; % 3x4 calibration matrix for right camera A*[R,T], where om is rotation vector and rodrigues() transform it into rotation matrix 
        Calib_P2 = struct('focus', fc1, 'center', cc1, 'distort', kc1, 'alpha_c', alpha_c1, 'matrix', C1);
        Calib_P1 = struct('focus', fc2, 'center', cc2, 'distort', kc2, 'alpha_c', alpha_c2, 'matrix', C2);
    end
    
    if strcmp(setName,'set1')
        FrmNum = 1634;
    elseif strcmp(setName,'set2')
        FrmNum = 1573;
    elseif strcmp(setName,'set3')
        FrmNum = 899;
    elseif strcmp(setName,'set4') 
        FrmNum = 2427;
    elseif strcmp(setName,'set5') 
        FrmNum = 3389; 
    end
    
function [out,dout]=rodrigues(in)
[m,n] = size(in);
bigeps = 10e+20*eps;
if ((m==1)&&(n==3))||((m==3)&&(n==1))%% it is a rotation vector
    theta = norm(in);
    if theta < eps
        R = eye(3);
        dRdin = [0 0 0;
            0 0 1;
            0 -1 0;
            0 0 -1;
            0 0 0;
            1 0 0;
            0 1 0;
            -1 0 0;
            0 0 0];
    else
        if n==length(in); in=in'; end; 	%% make it a column vec. if necess.
        dm3din = [eye(3);in'/theta];
        omega = in/theta;
        dm2dm3 = [eye(3)/theta -in/theta^2; zeros(1,3) 1];
        alpha = cos(theta);
        beta = sin(theta);
        gamma = 1-cos(theta);
        omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        A = omega*omega';
        
        dm1dm2 = zeros(21,4);
        dm1dm2(1,4) = -sin(theta);
        dm1dm2(2,4) = cos(theta);
        dm1dm2(3,4) = sin(theta);
        dm1dm2(4:12,1:3) = [0 0 0 0 0 1 0 -1 0;
            0 0 -1 0 0 0 1 0 0;
            0 1 0 -1 0 0 0 0 0]';

        w1 = omega(1);
        w2 = omega(2);
        w3 = omega(3);

        dm1dm2(13:21,1) = [2*w1;w2;w3;w2;0;0;w3;0;0];
        dm1dm2(13: 21,2) = [0;w1;0;w1;2*w2;w3;0;w3;0];
        dm1dm2(13:21,3) = [0;0;w1;0;0;w2;w1;w2;2*w3];

        R = eye(3)*alpha + omegav*beta + A*gamma;

        dRdm1 = zeros(9,21);
        dRdm1([1 5 9],1) = ones(3,1);
        dRdm1(:,2) = omegav(:);
        dRdm1(:,4:12) = beta*eye(9);
        dRdm1(:,3) = A(:);
        dRdm1(:,13:21) = gamma*eye(9);

        dRdin = dRdm1 * dm1dm2 * dm2dm3 * dm3din;
    end;
    out = R;
    dout = dRdin;
    %% it is prob. a rot matr.
elseif ((m==n)&&(m==3)&&(norm(in'*in-eye(3))<bigeps)&&(abs(det(in)-1)< bigeps))
    R = in;
    % project the rotation matrix to SO(3);
    [U,S,V] = svd(R);
    R = U*V';
    tr = (trace(R)-1)/2;
    dtrdR = [1 0 0 0 1 0 0 0 1]/2;
    theta = real(acos(tr));
    if sin(theta) >= 1e-4,
        dthetadtr = -1/sqrt(1-tr^2);
        dthetadR = dthetadtr * dtrdR;
        vth = 1/(2*sin(theta));
        dvthdtheta = -vth*cos(theta)/sin(theta);
        dvar1dtheta = [dvthdtheta;1];
        dvar1dR =  dvar1dtheta * dthetadR;
        om1 = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
        dom1dR = [0 0 0 0 0 1 0 -1 0;
            0 0 -1 0 0 0 1 0 0;
            0 1 0 -1 0 0 0 0 0];
        % var = [om1;vth;theta];
        dvardR = [dom1dR;dvar1dR];
        % var2 = [om;theta];
        om = vth*om1;
        domdvar = [vth*eye(3) om1 zeros(3,1)];
        dthetadvar = [0 0 0 0 1];
        dvar2dvar = [domdvar;dthetadvar];
        out = om*theta;
        domegadvar2 = [theta*eye(3) om];
        dout = domegadvar2 * dvar2dvar * dvardR;
    else
        if tr > 0; 			% case norm(om)=0;
            out = [0 0 0]';
            dout = [0 0 0 0 0 1/2 0 -1/2 0;
                0 0 -1/2 0 0 0 1/2 0 0;
                0 1/2 0 -1/2 0 0 0 0 0];
        else
            % case norm(om)=pi;
            if(0)
                %% fixed April 6th by Bouguet -- not working in all cases!
                out = theta * (sqrt((diag(R)+1)/2).*[1;2*(R(1,2:3)>=0)'-1]);
                %keyboard;
            else
                % Solution by Mike Burl on Feb 27, 2007
                % This is a better way to determine the signs of the
                % entries of the rotation vector using a hash table on all
                % the combinations of signs of a pairs of products (in the
                % rotation matrix)
                % Define hashvec and Smat
                hashvec = [0; -1; -3; -9; 9; 3; 1; 13; 5; -7; -11];
                Smat = [1,1,1; 1,0,-1; 0,1,-1; 1,-1,0; 1,1,0; 0,1,1; 1,0,1; 1,1,1; 1,1,-1;
                    1,-1,-1; 1,-1,1];
                M = (R+eye(3,3))/2;
                uabs = sqrt(M(1,1));
                vabs = sqrt(M(2,2));
                wabs = sqrt(M(3,3));
                mvec = [M(1,2), M(2,3), M(1,3)];
                syn  = ((mvec > 1e-4) - (mvec < -1e-4)); % robust sign() function
                hash = syn * [9; 3; 1];
                idx = hash == hashvec;
                svec = Smat(idx,:)';
                out = theta * [uabs; vabs; wabs] .* svec;

            end;

            if nargout > 1,
                fprintf(1,'WARNING!!!! Jacobian domdR undefined!!!\n');
                dout = NaN*ones(3,9);
            end;
        end;
    end;
else
    error('Neither a rotation matrix nor a rotation vector were provided');
end;