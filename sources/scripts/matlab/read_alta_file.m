

function [alta_file] = read_alta_file(filename)
%
% filename  : absolute path to an ALTA file
% alta_file : A structure that represents the DATA + Headers of the alta
% file
%
% alta_file structure :
%   alta_file.header : information from the header file
%   alta_file.data   : the numeric values of the data i.e., the matrix
%
% TODO: find a way to factor the common code regarding the Format of ALTA
% and the rest


alta_file = struct;


[header_ok, nbl_header,alta_file] = read_alta_file_header( alta_file, filename );

if header_ok
    [alta_file,data_ok] = read_alta_file_data( alta_file, filename, nbl_header );
end

if ~header_ok || ~data_ok
    disp('COULD NOT RETRIEVE CORRECT OR DATA HEADER...aborting...');
    return;
end


end

function [alta_file, reading_ok]=read_alta_file_data(alta_file, filename, nbl_header)

    imported_data_struct = importdata(filename,' ',nbl_header);
    
    %TODO Check that the size of the matrix
    [nbl,nbc] = size(imported_data_struct.data);
    if nbc == sum(alta_file.header.dim)
        alta_file.data = imported_data_struct.data;
        reading_ok = true;
    else
        disp('Inconsistent dimensions between headers and data');
        reading_ok = false;
    end

end


function [header_ok, nbl_header,alta_file] = read_alta_file_header(alta_file, filename )
%
%
% header_ok true if the header was correctly read
% nbl_header : number of lines detected in the header
%

[fid,errmsg] = fopen(filename);
if fid == -1
    disp(errmsg);
    header_ok = false;
    return;
end

input_param_map  =  alta_inputParametrizationNameMap();
output_param_map = alta_outputParametrizationNameMap();


nbl_header = 0;
tline = fgetl(fid);
while strfind(tline, '#')
    
    line_tokens = textscan(tline,'%s');
    
    alta_file = adjustHeaderInformationFromLineTokens( alta_file, line_tokens);
    
    tline=fgetl(fid);
    nbl_header = nbl_header + 1;
end

fclose(fid);
header_ok = true;

end

function [alta_file,parsing_ok]=adjustHeaderInformationFromLineTokens( alta_file, line_tokens)
    
parsing_ok = true;
first_token = cellstr(line_tokens{1}(1));
    
[alta_data_header_start, alta_data_header_stop]= alta_dataHeaderTokens();
if strcmp( first_token, alta_data_header_start) || strcmp(first_token,alta_data_header_stop)
    return;
end


if strcmp(first_token, dimensionToken())
    dim_in  = str2double( line_tokens{1}(2) );
    dim_out = str2double( line_tokens{1}(3) );   
    alta_file.header.dim = [ dim_in dim_out];
    return;
end

if strcmp(first_token, alta_vsToken() )
    
    [nb_tokens,~] = size(line_tokens{1});
    vs_numbers = zeros(nb_tokens-1,1);
    
    for i=2:nb_tokens
        vs_numbers(i-1,1) = str2double( line_tokens{1}(i) );
    end
    alta_file.header.VS = vs_numbers;
    return;
end

input_param_map = alta_inputParametrizationNameMap();

[param_in_token,param_out_token] = alta_parametrizationTokens();

if strcmp(first_token, param_in_token )
    
    param_in_name = line_tokens{1}(2);
    
    if isKey(input_param_map,param_in_name)
        alta_file.header.param_in = param_in_name;
    else
        parsing_ok = false;
    end
end

output_param_map = alta_outputParametrizationNameMap();

if strcmp(first_token, param_out_token )
    param_out_name = line_tokens{1}(2);
    if isKey(output_param_map,param_out_name)
        alta_file.header.param_out = param_out_name;
    else
        parsing_ok = false;
    end
 
end


%disp('We had a comment...skipping line');


end

function   input_param_map = alta_inputParametrizationNameMap()
input_param_keys =  {'RUSIN_TH_PH_TD_PD',...
                    'RUSIN_TH_PH_TD',...
                    'RUSIN_TH_TD_PD', ...
                    'RUSIN_TH_TD', ...   
                    'RUSIN_VH_VD', ...  
                    'RUSIN_VH', ...      
                    'COS_TH_TD', ...            
                    'COS_TH', ...
                    'SCHLICK_TK_PK', ...       
                    'SCHLICK_VK', ...          
                    'SCHLICK_TL_TK_PROJ_DPHI', ...
                    'COS_TK', ...                
                    'RETRO_TL_TVL_PROJ_DPHI', ...
                    'STEREOGRAPHIC', ...        
                    'SPHERICAL_TL_PL_TV_PV', ...
                    'COS_TLV', ...               
                    'COS_TLR', ...               
                    'ISOTROPIC_TV_TL', ...       
                    'ISOTROPIC_TV_TL_DPHI', ...                                
                    'ISOTROPIC_TV_PROJ_DPHI',...                             
                    'ISOTROPIC_TL_TV_PROJ_DPHI',...
                    'ISOTROPIC_TD_PD',...      
                    'BARYCENTRIC_ALPHA_SIGMA',...
                    'CARTESIAN',...   
                    'UNKNOWN_INPUT'};
       
input_param_values = ones(size(input_param_keys));
input_param_map = containers.Map(input_param_keys, input_param_values );             
end


function  output_param_map = alta_outputParametrizationNameMap()
output_param_keys = {'INV_STERADIAN',...              
                     'INV_STERADIAN_COSINE_FACTOR',...                                                    
                     'ENERGY',...
                     'RGB_COLOR',...
                     'XYZ_COLOR',...
                     'UNKNOWN_OUTPUT'};
[~,nbe] = size(output_param_keys);
output_param_values = ones(1, nbe);
output_param_map  = containers.Map(output_param_keys, output_param_values );

end


function   [alta_data_header_start, alta_data_header_stop]= alta_dataHeaderTokens()
alta_data_header_start = '#ALTA_HEADER_DATA';
alta_data_header_stop  = '#ALTA_HEADER_END'; 
end


function  dim_token=dimensionToken()
dim_token = '#DIM'; 
end

function [param_in_token,param_out_token]=alta_parametrizationTokens()
param_in_token  = '#PARAM_IN';
param_out_token = '#PARAM_OUT';
end

function vs_token=alta_vsToken()
vs_token ='#VS';
end

function binary_data_token = alta_binaryToken()
binary_data_token = '#BINARY';
end


