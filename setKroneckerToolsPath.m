%  The KroneckerTools can be obtained from
%       https://github.com/jborggaard/KroneckerTools
%
%  Download KroneckerTools, then adjust the path variable below if needed.
KroneckerToolsPath = '/Volumes/borggaard4/Software/MyPublicSoftware/KroneckerTools';

addpath([KroneckerToolsPath,'/src'])
addpath([KroneckerToolsPath,'/util'])

if ( exist([KroneckerToolsPath,'/tensor_recursive'],'dir') )
  addpath([KroneckerToolsPath,'/tensor_recursive'])
end
