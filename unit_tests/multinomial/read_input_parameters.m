function inputStructure = read_input_parameters(filename) 

fid = fopen(filename,'r');
data = textscan(fid,'%s%f','Delimiter','=');
fclose(fid);
for kk = 1:length(data{1})
    inputStructure.(strtrim(data{1}{kk})) = data{2}(kk);
end
