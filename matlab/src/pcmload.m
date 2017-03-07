function PCM = pcmload(NCFILE)
% pcmload Load a Profile Classification Model (PCM) from a netcdf file
%
% PCM = pcmload(NCFILE)
%
% See Also: pcmsave

% PCM is Profile Classification Modelling
% Copyright (C) 2016-2017, OBIDAM Developpers
% For more information, see http://framagit.org/obidam/pcm
% Created: 2016-04-14 (G. Maze, Ifremer, Laboratoire d'Océanographie Physique et Spatiale)
% Revised: 2016-12-09 (G. Maze) Added netcdf functions xtype and listAtt

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

[PATHSTR,NAME,EXT] = fileparts(NCFILE);
if isempty(EXT)
	NCFILE = sprintf('%s.nc',NCFILE);
end% if 

% Open Netcdf file:
scope = netcdf.open(NCFILE,'netcdf4');
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');

% Get PCM File format version:
try
	ver = netcdf.getAtt(scope,NC_GLOBAL,'format_version');
catch
	try 
		ver = netcdf.getAtt(scope,NC_GLOBAL,'FormatVersion');
	catch
		error('I dont understand the Format Version parameter of this file !');
	end
end

% Load data according to format version
try
	switch ver
		case '1.0'
			PCM = loadthis10(scope);
		otherwise
			error(sprintf('I don''t how to load this PCM file format: %s',ver));
	end% switch 
	PCM = orderfields(PCM);	
catch ME
	netcdf.close(scope)
	rethrow(ME)
end%try

end %functionpcmload

%- Load file format version 1.0
function PCM = loadthis10(scope)
	
	getvar = @(sc, vn, vt) netcdf.getVar(sc,netcdf.inqVarID(sc,vn),vt);
	
	PCM = struct();
	
	% Define useful constants:
	fillValue = -9999;
	NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');	
	
	%-- Load global attributes
	[GbAttnames GbAttIDs GbAttNCType] = nclistAtt(scope,'NC_GLOBAL');
	
	% Load global attributes variables:
	plist = {'doREDUCE','normalization','Np','readme'};
	for ip = 1 : length(plist)
		try
			[~,iat] = intersect(GbAttnames,sprintf('PCM_%s',plist{ip}));
			mtype = ncxtype(ncxtype(GbAttNCType{iat}),'matlab');
			att = netcdf.getAtt(scope,NC_GLOBAL,GbAttnames{iat},mtype);
			PCM = setfield(PCM,plist{ip},att);
		catch ME
			rethrow(ME)
		end	
	end% for ip	
	PCM.doREDUCE  = eval(PCM.doREDUCE);
	
	%-- Load model vertical dimensions / axis:
	PCM.DPTmodel = netcdf.getVar(scope,netcdf.inqVarID(scope,'DEPTH_MODEL'),'double');
	Nz = length(PCM.DPTmodel);
	
	%-- Load Normalization data:
	norm_scope = netcdf.inqNcid(scope,'Normalization');		
	switch PCM.normalization
		case 0 %--- No norm: Xn = X;
			% Nothing to add the netcdf file
			
		case 1 %--- Center/standardize: Xn = (X - repmat(X_ave,[1 Np]))./repmat(X_std,[1 Np]);
			% Sample mean of profiles:
			PCM.X_ave = getvar(norm_scope,'X_ave','double');
			% Sample std of profiles:
			PCM.X_std = getvar(norm_scope,'X_std','double');
			
		case 2 %--- Center only: Xn = (X - repmat(X_ave,[1 Np]));
			% Sample mean of profiles:
			PCM.X_ave = getvar(norm_scope,'X_ave','double');
		
	end% switch 
	
	%-- Load Reduction data (PCA based)
	reduction_scope = netcdf.inqNcid(scope,'Reduction');
	
	switch PCM.doREDUCE
		case 0 %--- No Reduction
			
		case 1 %--- Reduction with PCA
			PCM.X_ref = getvar(reduction_scope,'X_ref','double');
			PCM.EOFs = getvar(reduction_scope,'EOFs','double');
			PCM.V = getvar(reduction_scope,'EOFvar','double')';			
			PCM.maxvar = sum(PCM.V)*100;
	end% switch 

	%-- Load GMM variables
	cmodel_scope = netcdf.inqNcid(scope,'ClassificationModel');
	[~, PCM.K] = netcdf.inqDim(cmodel_scope,netcdf.inqDimID(cmodel_scope,'CLASS'));
	PCM.covarTYPE = netcdf.getAtt(cmodel_scope,NC_GLOBAL,'Covariance_matrix_original_form');

	if PCM.doREDUCE
		Nd = size(PCM.EOFs,2);	
	else
		Nd = size(PCM.DPTmodel,1);
	end% if
	
	PCM.mix = struct(); % Create a Netlab GMM mixture structure
	PCM.mix.type = 'gmm';
	PCM.mix.nin = Nd;
	PCM.mix.ncentres = PCM.K;
	PCM.mix.covar_type = PCM.covarTYPE;
		
	PCM.mix.priors  = getvar(cmodel_scope,'priors','double')';
	PCM.mix.centres = getvar(cmodel_scope,'centers','double')';
	PCM.mix.covars  = getvar(cmodel_scope,'covariances','double');
	PCM.LLH  = getvar(cmodel_scope,'llh','double')';

	% Update Covariance matrix according to covar_type and Netlab way to store it:
	% These info can be found into netlab gmm.m function
	% rq: nwts is the Nb of GMM unknowns
	switch PCM.mix.covar_type
		case 'full'
			covars = PCM.mix.covars;
			PCM.mix.nwts = PCM.mix.ncentres + PCM.mix.ncentres*PCM.mix.nin + PCM.mix.ncentres*PCM.mix.nin*PCM.mix.nin;
		case 'diag'
			covars = zeros(PCM.K, Nd);
			for ik = 1 : PCM.K
				covars(ik,:) = diag(PCM.mix.covars(:,:,ik));
			end% for ik
			PCM.mix.nwts = PCM.mix.ncentres + PCM.mix.ncentres*PCM.mix.nin + PCM.mix.ncentres*PCM.mix.nin;
		case 'spherical'
			covars = ones(1, PCM.K);
			for ik = 1 : PCM.K
				PCM.mix.covars(:,:,ik)
				covars(1,ik) = PCM.mix.covars(1,1,ik);
			end% for ik
			PCM.mix.nwts = PCM.mix.ncentres + PCM.mix.ncentres*PCM.mix.nin + PCM.mix.ncentres;		  
	end% switch	
	PCM.mix.covars = covars;
	
end%end function

function varargout = nclistAtt(ncid,varid)
% listAtt List all Attributes of a variable
%
% [AttName, AttIDs, AttNCtype] = listAtt(ncid,var[id])
% 
% List all Attributes of a variable. If no output required displau list on screen
%
% Inputs:
%       ncid: netcdf scope of the variable
%       var[id]: Variable name (string) or ID (int). 
%               It can also be 'NC_GLOBAL' for global variables
%
% Outputs:
%       AttName: Cell of Attribute's name (string)
%       AttIDs: Attribute's IDs (int)
%       AttXtype: Cell of Attribute's NC_TYPE (string)
% 
% Eg:
%       listAtt(ncid,3)
%       listAtt(ncid,'TEMP')
%       listAtt(ncid,'GLOBAL')
% 
% See Also: listVar
%
% Copyright (c) 2016, Guillaume Maze (Ifremer, Laboratoire d'Océanographie Physique et Spatiale).
% For more information, see the http://codes.guillaumemaze.org
% Created: 2016-04-14 (G. Maze)

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%       * Redistributions of source code must retain the above copyright notice, this list of 
%       conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright notice, this list 
%       of conditions and the following disclaimer in the documentation and/or other materials 
%       provided with the distribution.
%       * Neither the name of the Ifremer, Laboratoire d'Océanographie Physique et Spatiale nor the names of its contributors may be used 
%       to endorse or promote products derived from this software without specific prior 
%       written permission.
%
% THIS SOFTWARE IS PROVIDED BY Guillaume Maze ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, 
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Guillaume Maze BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if ischar(varid)
        if strcmp(varid,'NC_GLOBAL') | strcmp(varid,'GLOBAL')
                varid = netcdf.getConstant('NC_GLOBAL');
        else
                varid = netcdf.findVarID(ncid,varid);
        end% if 
end% if

if ischar(ncid)
        error('Invalid ncid')
end% if 

trashed = false;
attid = 0;%keyboard
while ~trashed
        try
                attname = netcdf.inqAttName(ncid,varid,attid);
                [xtype,attlen] = netcdf.inqAtt(ncid,varid,attname);
                attype = ncxtype(xtype);
                it = attid+1;
                Anames{it} = attname;
                Aids{it} = attid;
                Atype{it} = attype;
                dstr = sprintf('\t#%3.1d: %20s [%s]',attid,attname,attype);
                RESdisp(it) = {dstr};
        catch   
                trashed = true;
        end% try
        attid = attid + 1;
end% for it

switch nargout
        case 1
                varargout(1) = {Anames};
        case 2
                varargout(1) = {Anames};
                varargout(2) = {Aids};
        case 3
                varargout(1) = {Anames};
                varargout(2) = {Aids};
                varargout(3) = {Atype};
        otherwise
                if varid == netcdf.getConstant('NC_GLOBAL')
                        varname = 'NC_GLOBAL';
                        disp(sep('-',sprintf(' LIST OF GLOBAL ATTRIBUTE(S) ')))                 
                else
                        varname = netcdf.inqVar(ncid,varid);
                        disp(sep('-',sprintf(' LIST OF ATTRIBUTE(S) FOR %s ',varname)))                 
                end% if 
                disp(sprintf('\t#IDS: %20s [%s]','ATTRIBUTE''S NAME','TYPE'))
                sep
                for iv = 1 : it
                        disp(RESdisp{iv});
                end% for iv
                sep
end% switch 


end %functionlistAtt

function varargout = ncxtype(xt,varargin)
% xtype Corresponding value between NC_TYPE and XTYPE
%
% NC_TYPE = xtype(XTYPE) Return the NC_TYPE corresponding to XTYPE
% XTYPE = xtype(NC_TYPE) Return the XTYPE corresponding to NC_TYPE
%
% MATLAB_TYPE = xtype(XTYPE,'matlab') Return the MATLAB_TYPE corresponding to XTYPE
% XTYPE = xtype(MATLAB_TYPE) Return the XTYPE corresponding to MATLAB_TYPE
%
% XTYPE are the set of predefined netCDF external data types. It can be:
% NC_BYTE, NC_CHAR, NC_SHORT, NC_INT, NC_FLOAT and NC_DOUBLE
%
% Copyright (c) 2016, Guillaume Maze (Ifremer, Laboratoire d'Océanographie Physique et Spatiale).
% For more information, see the http://codes.guillaumemaze.org
% Created: 2016-04-14 (G. Maze)

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%       * Redistributions of source code must retain the above copyright notice, this list of 
%       conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright notice, this list 
%       of conditions and the following disclaimer in the documentation and/or other materials 
%       provided with the distribution.
%       * Neither the name of the Ifremer, Laboratoire d'Océanographie Physique et Spatiale nor the names of its contributors may be used 
%       to endorse or promote products derived from this software without specific prior 
%       written permission.
%
% THIS SOFTWARE IS PROVIDED BY Guillaume Maze ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, 
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Guillaume Maze BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

NC_TYPES  = {'NC_BYTE','NC_CHAR','NC_SHORT','NC_INT','NC_FLOAT','NC_DOUBLE'};
MAT_TYPES = {'int8',   'char',   'int16',   'int32', 'single',  'double'};
for il = 1 : length(NC_TYPES)
        XTYPE(il) = netcdf.getConstant(NC_TYPES{il});
end% for lt 

if ischar(xt) % NC_TYPES -> XTYPE
        if nargin==2 & strcmp(lower(varargin{1}),'matlab')
                [~,ix] = intersect(MAT_TYPES,xt);
                if isempty(ix)
                        error('Unknown Matlab TYPE (try with lower case or Netcdf type)')
                end% if
        else    
                [~,ix] = intersect(NC_TYPES,xt);
                if isempty(ix)
                        error('Unknown Netcdf TYPE (try with upper case or Matlab type)')
                end% if 
        end% if 
        varargout(1) = {XTYPE(ix)};
else % XTYPE -> NC_TYPES
        varargout(1) = {NC_TYPES{xt}};
        if nargin==2 & strcmp(lower(varargin{1}),'matlab')
                varargout(1) = {MAT_TYPES{xt}};
        end% if 
end% if 

end %functionxtype