function varargout = savePCM(NCFILE,PCM)
% savePCM Save a Profile Classification Model (PCM) into a netcdf file
%
% [] = savePCM(NCFILE,PCM)
%
% Save PCM Profile Classification Model structure into a netcdf file NCFILE
%
% See Also: pcmload
%
% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://github.com/obidam/pcm
% Created: 2016-04-13 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)

% This file is part of OBIDAM/PCM.
%     OBIDAM/PCM is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%     OBIDAM/PCM is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%     You should have received a copy of the GNU General Public License
%     along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

% File format version
ver = '1.0';

[PATHSTR,NAME,EXT] = fileparts(NCFILE);
if isempty(EXT)
	NCFILE = sprintf('%s.nc',NCFILE);
end% if 

% Check if file exist
reply = 'Y';
if exist(NCFILE,'file')
	reply = input('This model file already exists, do you want to overwite it? Y/N [Y]:','s');
	if strcmp(lower(reply),'n'), 
		reply = 'N'; 
	else
		reply = 'Y'; % Default response
	end% if 
end% if
if strcmp(reply,'Y')
	% Delete file to proceed
	delete(NCFILE);
else
	% Abort saving
	error(sprintf('The PCM was not saved in %s because it already exists',NCFILE))
end% if

% Proceed to saving
switch ver
	case '1.0' %- Version 1.0
		
		%-- Read and check usefull sizes from the PCM structure
		Nz = length(PCM.DPTmodel);
		
		if PCM.doREDUCE
			if size(PCM.EOFs,2) ~= PCM.mix.nin
				error('I can''t save this PCM because the PCA eigenvectors size doesn''t match the GMM dimensions');
			else
				Nd = size(PCM.EOFs,2);	
			end% if 
		else
			Nd = Nz;
		end% if 
		
		if PCM.K ~= PCM.mix.ncentres
			error('I can''t save this PCM because the number of class doesn''t match the GMM structure');
		else
			K = PCM.K;
		end% if 
		
		% Open netcdf4 file:
		scope = netcdf.create(NCFILE,'netcdf4');
		
		% Define useful stuff:
		NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
		fillValue = -9999;
		dimidVAL = netcdf.defDim(scope,'ONEVAL',1);
		
		%-- Insert global attributes
		netcdf.putAtt(scope,NC_GLOBAL,'title','PCM parameters')
		netcdf.putAtt(scope,NC_GLOBAL,'long_title','Profile Classification Model parameters for prediction')
		netcdf.putAtt(scope,NC_GLOBAL,'institution','Ifremer/LOPS')
		
		netcdf.putAtt(scope,NC_GLOBAL,'software','Profile Classification Model - Matlab Toolbox (c) Ifremer')
		netcdf.putAtt(scope,NC_GLOBAL,'software_version','1.0')
		netcdf.putAtt(scope,NC_GLOBAL,'software_source','')
		netcdf.putAtt(scope,NC_GLOBAL,'format_version',ver)
		
		bool = {'false','true'};
		netcdf.putAtt(scope,NC_GLOBAL,'PCM_doREDUCE',bool{PCM.doREDUCE+1});
		netcdf.putAtt(scope,NC_GLOBAL,'PCM_normalization',int32(PCM.normalization));
		%netcdf.putAtt(scope,NC_GLOBAL,'PCM_K',int32(PCM.K));
		netcdf.putAtt(scope,NC_GLOBAL,'PCM_Np',int32(PCM.Np));
		%netcdf.putAtt(scope,NC_GLOBAL,'PCM_covarTYPE',PCM.covarTYPE);
		netcdf.putAtt(scope,NC_GLOBAL,'PCM_readme',PCM.readme);
		
		netcdf.putAtt(scope,NC_GLOBAL,'Conventions','CF-1.6')   
		netcdf.putAtt(scope,NC_GLOBAL,'Conventions_help','http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html')

		netcdf.putAtt(scope,NC_GLOBAL,'CreationDate',datestr(now,'yyyy/mm/dd HH:MM:SS'))
		netcdf.putAtt(scope,NC_GLOBAL,'CreatedBy',getenv('LOGNAME'))

		%-- Save variables required for Interpolation
		%interp_scope = netcdf.defGrp(scope,'Interpolation');
				
		% Define model vertical dimensions / axis:
		dimidZ = netcdf.defDim(scope,'DEPTH_MODEL',Nz);
		varid  = netcdf.defVar(scope,'DEPTH_MODEL','double',[dimidZ]);
		netcdf.putAtt(scope,varid,'standard_name','depth');
		netcdf.putAtt(scope,varid,'long_name','Vertical distance below the surface');
		netcdf.putAtt(scope,varid,'units','m');
		netcdf.putAtt(scope,varid,'comment','Depth axis used by the classification model');
		netcdf.defVarFill(scope,varid,false,fillValue);
		netcdf.putVar(scope,varid,PCM.DPTmodel);
		     			
		%-- Save variables required for Normalization
		norm_scope = netcdf.defGrp(scope,'Normalization');
		
		switch PCM.normalization
			case 0 %--- No norm: Xn = X;
				% Nothing to add the netcdf file
				netcdf.putAtt(norm_scope,NC_GLOBAL,'Comment','The model doesn''t not use normalization')
				netcdf.putAtt(norm_scope,NC_GLOBAL,'normalization',int32(PCM.normalization))
				
			case 1 %--- Center/standardize: Xn = (X - repmat(X_ave,[1 Np]))./repmat(X_std,[1 Np]);
				netcdf.putAtt(norm_scope,NC_GLOBAL,'Comment','The model centers and standardizes data at each depth levels')
				netcdf.putAtt(norm_scope,NC_GLOBAL,'normalization',int32(PCM.normalization))
				
				% Sample mean of profiles:
				varid = netcdf.defVar(norm_scope,'X_ave','double',[dimidZ]);
				netcdf.putAtt(norm_scope,varid,'name','X_ave');
				netcdf.putAtt(norm_scope,varid,'long_name','Sample mean of profiles');
				netcdf.putAtt(norm_scope,varid,'units','');   
				%netcdf.putAtt(norm_scope,varid,'comment',sprintf('Used for normalization=%i',PCM.normalization));
				netcdf.defVarFill(norm_scope,varid,false,fillValue);
			    netcdf.putVar(norm_scope,varid,PCM.X_ave);

				% Sample std of profiles:
				varid = netcdf.defVar(norm_scope,'X_std','double',[dimidZ]);
				netcdf.putAtt(norm_scope,varid,'name','X_std');
				netcdf.putAtt(norm_scope,varid,'long_name','Sample standard deviation of profiles');
				netcdf.putAtt(norm_scope,varid,'units','');
				%netcdf.putAtt(norm_scope,varid,'comment',sprintf('Used for normalization=%i',PCM.normalization));
				netcdf.defVarFill(norm_scope,varid,false,fillValue);
			    netcdf.putVar(norm_scope,varid,PCM.X_std);
				
			case 2 %--- Center only: Xn = (X - repmat(X_ave,[1 Np]));
				netcdf.putAtt(norm_scope,NC_GLOBAL,'Comment','The model centers data at each depth levels')
				netcdf.putAtt(norm_scope,NC_GLOBAL,'normalization',int32(PCM.normalization))
					
				% Sample mean of profiles:
				varid = netcdf.defVar(norm_scope,'X_ave','double',[dimidZ]);
				netcdf.putAtt(norm_scope,varid,'name','X_ave');
				netcdf.putAtt(norm_scope,varid,'long_name','Sample mean of profiles');
				netcdf.putAtt(norm_scope,varid,'units','');   
				%netcdf.putAtt(norm_scope,varid,'comment',sprintf('Used for normalization=%i',PCM.normalization));
				netcdf.defVarFill(norm_scope,varid,false,fillValue);
			    netcdf.putVar(norm_scope,varid,PCM.X_ave);
			
		end% switch 
					
		%-- Save variables required for compression (PCA based)
		reduction_scope = netcdf.defGrp(scope,'Reduction');
		
		switch PCM.doREDUCE
			case 0 %--- No Reduction
				netcdf.putAtt(reduction_scope,NC_GLOBAL,'Comment','The model does not use dimensionality reduction')
				dimidD = dimidZ;
				
			case 1 %--- Reduction with PCA
				netcdf.putAtt(reduction_scope,NC_GLOBAL,'Comment','The model uses dimensionality reduction with PCA')

				% Define netcdf reduced dimensions:
				dimidD = netcdf.defDim(scope,'REDUCED_DIM',Nd);
						
				% Reference profile for reduction:
				varid = netcdf.defVar(reduction_scope,'X_ref','double',[dimidZ]);
				netcdf.putAtt(reduction_scope,varid,'name','X_ref');
				netcdf.putAtt(reduction_scope,varid,'long_name','Reference profile for PCA reduction');
				netcdf.putAtt(reduction_scope,varid,'comment','Note that if normalization=1/2, data are already centered before reduction so X_ref is epsilon');
				netcdf.putAtt(reduction_scope,varid,'units','');
				netcdf.defVarFill(reduction_scope,varid,false,fillValue);
			    netcdf.putVar(reduction_scope,varid,PCM.X_ref);

				% Define netcdf reduced axis:
				varid = netcdf.defVar(reduction_scope,'EOFs','double',[dimidZ dimidD]);         
				netcdf.putAtt(reduction_scope,varid,'name','Reduced space');
				netcdf.putAtt(reduction_scope,varid,'long_name','Eigenvectors for PCA dimensionality reduction');
				netcdf.putAtt(reduction_scope,varid,'howto','X_reduc = EOFs''*(X-X_ref);')
				netcdf.putAtt(reduction_scope,varid,'units','');
				netcdf.defVarFill(reduction_scope,varid,false,fillValue);
				netcdf.putVar(reduction_scope,varid,PCM.EOFs);

				% Define explained variance for each reduced axis:
				varid = netcdf.defVar(reduction_scope,'EOFvar','double',[dimidD]);         
				netcdf.putAtt(reduction_scope,varid,'name','Percentage of variance');
				netcdf.putAtt(reduction_scope,varid,'long_name','Percentage of variance attributed to each PCA dimension');
				netcdf.putAtt(reduction_scope,varid,'units','%');
				netcdf.defVarFill(reduction_scope,varid,false,fillValue);
				netcdf.putVar(reduction_scope,varid,PCM.V);
				
		end% switch 
				
		%-- Save GMM variables
		cmodel_scope = netcdf.defGrp(scope,'ClassificationModel');
		netcdf.putAtt(cmodel_scope,NC_GLOBAL,'Comment','Parameters of the Gaussian Mixture Model used by the PCM')
		%netcdf.putAtt(cmodel_scope,NC_GLOBAL,'Number_of_Class',int32(PCM.K));
		netcdf.putAtt(cmodel_scope,NC_GLOBAL,'Covariance_matrix_original_form',PCM.mix.covar_type)				
				
		% Define netcdf GMM class "dimensions":
		dimidK = netcdf.defDim(cmodel_scope,'CLASS',K);		

		% Log likelihood
		varid = netcdf.defVar(cmodel_scope,'llh','double',[dimidVAL]);         
		netcdf.putAtt(cmodel_scope,varid,'name','llh');
		netcdf.putAtt(cmodel_scope,varid,'long_name','Negative Log likelihood');    
		netcdf.putAtt(cmodel_scope,varid,'howto','-sum(log(gmmprob(mix,x)))');
		netcdf.putAtt(cmodel_scope,varid,'comment','The smaller the better');
		netcdf.putAtt(cmodel_scope,varid,'units','');
		netcdf.defVarFill(cmodel_scope,varid,false,fillValue);
		netcdf.putVar(cmodel_scope,varid,PCM.LLH);
		
		% Priors
		varid = netcdf.defVar(cmodel_scope,'priors','double',[dimidK]);         
		netcdf.putAtt(cmodel_scope,varid,'name','Priors');
		netcdf.putAtt(cmodel_scope,varid,'long_name','Class weights');    
		netcdf.putAtt(cmodel_scope,varid,'comments','sum(priors) = 1 & 0 < priors < 1');
		netcdf.putAtt(cmodel_scope,varid,'units','');
		netcdf.defVarFill(cmodel_scope,varid,false,fillValue);
		netcdf.putVar(cmodel_scope,varid,PCM.mix.priors);
		
		% Gaussian Centers
		varid = netcdf.defVar(cmodel_scope,'centers','double',[dimidD dimidK]);         
		netcdf.putAtt(cmodel_scope,varid,'name','Centers');
		netcdf.putAtt(cmodel_scope,varid,'long_name','Gaussian Class centers');
		netcdf.putAtt(cmodel_scope,varid,'units','');
		netcdf.defVarFill(cmodel_scope,varid,false,fillValue);
		netcdf.putVar(cmodel_scope,varid,PCM.mix.centres');
		
		% Gaussian covariances
		switch PCM.mix.covar_type
			case 'full'
				cf = 'plain';
				covars = PCM.mix.covars;
			case 'diag'	
				cf = 'diagonal';
				% Recompose a squared covariance matrix:
				covars = zeros(Nd, Nd, K);
				for ik = 1 : K
					covars(:,:,ik) = diag(PCM.mix.covars(ik,:));
				end% for ik		
			case 'spherical'
				cf = 'spherical';
				% Recompose a squared covariance matrix:
				covars = zeros(Nd, Nd, K);
				for ik = 1 : K
					covars(:,:,ik) = PCM.mix.covars(ik)*eye(Nd);
				end% for ik				
		end% switch	
		varid = netcdf.defVar(cmodel_scope,'covariances','double',[dimidD dimidD dimidK]);         
		netcdf.putAtt(cmodel_scope,varid,'name','Covariances');
		netcdf.putAtt(cmodel_scope,varid,'long_name','Gaussian Class squared covariance matrices');
		netcdf.putAtt(cmodel_scope,varid,'units','');
		netcdf.defVarFill(cmodel_scope,varid,false,fillValue);
		netcdf.putVar(cmodel_scope,varid,covars);
				
		
		% Write and close file
		netcdf.close(scope);		
	
end% switch 


end %functionsavePCM


















