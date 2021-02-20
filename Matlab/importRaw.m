function MazzioJon3 = importRaw(filename, startRow, endRow)
%IMPORTFILE1 Import numeric data from a text file as a matrix.
%   MAZZIOJON3 = IMPORTFILE1(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   MAZZIOJON3 = IMPORTFILE1(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   MazzioJon3 = importfile1('Mazzio_Jon_3.asc', 8, 1920);
%
%    See also TEXTSCAN.


%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 8;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
MazzioJon3 = table;
MazzioJon3.timesecs = cell2mat(raw(:, 1));
MazzioJon3.VarName2 = cell2mat(raw(:, 2));
MazzioJon3.VarName3 = cell2mat(raw(:, 3));
MazzioJon3.VarName4 = cell2mat(raw(:, 4));
MazzioJon3.VarName5 = cell2mat(raw(:, 5));
MazzioJon3.VarName6 = cell2mat(raw(:, 6));
MazzioJon3.VarName7 = cell2mat(raw(:, 7));
MazzioJon3.VarName8 = cell2mat(raw(:, 8));
MazzioJon3.VarName9 = cell2mat(raw(:, 9));
MazzioJon3.VarName10 = cell2mat(raw(:, 10));
MazzioJon3.VarName11 = cell2mat(raw(:, 11));
MazzioJon3.VarName12 = cell2mat(raw(:, 12));
MazzioJon3.VarName13 = cell2mat(raw(:, 13));
MazzioJon3.VarName14 = cell2mat(raw(:, 14));
MazzioJon3.VarName15 = cell2mat(raw(:, 15));
MazzioJon3.VarName16 = cell2mat(raw(:, 16));
MazzioJon3.VarName17 = cell2mat(raw(:, 17));
MazzioJon3.VarName18 = cell2mat(raw(:, 18));
MazzioJon3.VarName19 = cell2mat(raw(:, 19));
MazzioJon3.VarName20 = cell2mat(raw(:, 20));
MazzioJon3.VarName21 = cell2mat(raw(:, 21));
MazzioJon3.VarName22 = cell2mat(raw(:, 22));
MazzioJon3.VarName23 = cell2mat(raw(:, 23));
MazzioJon3.VarName24 = cell2mat(raw(:, 24));
MazzioJon3.VarName25 = cell2mat(raw(:, 25));
MazzioJon3.VarName26 = cell2mat(raw(:, 26));
MazzioJon3.VarName27 = cell2mat(raw(:, 27));
MazzioJon3.VarName28 = cell2mat(raw(:, 28));
MazzioJon3.VarName29 = cell2mat(raw(:, 29));
MazzioJon3.VarName30 = cell2mat(raw(:, 30));
MazzioJon3.VarName31 = cell2mat(raw(:, 31));
MazzioJon3.VarName32 = cell2mat(raw(:, 32));
MazzioJon3.VarName33 = cell2mat(raw(:, 33));
MazzioJon3.VarName34 = cell2mat(raw(:, 34));
MazzioJon3.VarName35 = cell2mat(raw(:, 35));
MazzioJon3.VarName36 = cell2mat(raw(:, 36));
MazzioJon3.VarName37 = cell2mat(raw(:, 37));
MazzioJon3.VarName38 = cell2mat(raw(:, 38));
MazzioJon3.VarName39 = cell2mat(raw(:, 39));
MazzioJon3.VarName40 = cell2mat(raw(:, 40));
MazzioJon3.VarName41 = cell2mat(raw(:, 41));
MazzioJon3.VarName42 = cell2mat(raw(:, 42));
MazzioJon3.VarName43 = cell2mat(raw(:, 43));
MazzioJon3.VarName44 = cell2mat(raw(:, 44));
MazzioJon3.VarName45 = cell2mat(raw(:, 45));
MazzioJon3.VarName46 = cell2mat(raw(:, 46));
MazzioJon3.VarName47 = cell2mat(raw(:, 47));
MazzioJon3.VarName48 = cell2mat(raw(:, 48));
MazzioJon3.VarName49 = cell2mat(raw(:, 49));
MazzioJon3.VarName50 = cell2mat(raw(:, 50));
MazzioJon3.VarName51 = cell2mat(raw(:, 51));
MazzioJon3.VarName52 = cell2mat(raw(:, 52));
MazzioJon3.VarName53 = cell2mat(raw(:, 53));
MazzioJon3.VarName54 = cell2mat(raw(:, 54));
MazzioJon3.VarName55 = cell2mat(raw(:, 55));
MazzioJon3.VarName56 = cell2mat(raw(:, 56));
MazzioJon3.VarName57 = cell2mat(raw(:, 57));
MazzioJon3.VarName58 = cell2mat(raw(:, 58));
MazzioJon3.VarName59 = cell2mat(raw(:, 59));
MazzioJon3.VarName60 = cell2mat(raw(:, 60));
MazzioJon3.VarName61 = cell2mat(raw(:, 61));
MazzioJon3.VarName62 = cell2mat(raw(:, 62));
MazzioJon3.VarName63 = cell2mat(raw(:, 63));
MazzioJon3.VarName64 = cell2mat(raw(:, 64));
MazzioJon3.VarName65 = cell2mat(raw(:, 65));
MazzioJon3.VarName66 = cell2mat(raw(:, 66));
MazzioJon3.VarName67 = cell2mat(raw(:, 67));
MazzioJon3.VarName68 = cell2mat(raw(:, 68));
MazzioJon3.VarName69 = cell2mat(raw(:, 69));
MazzioJon3.VarName70 = cell2mat(raw(:, 70));
MazzioJon3.VarName71 = cell2mat(raw(:, 71));
MazzioJon3.VarName72 = cell2mat(raw(:, 72));
MazzioJon3.VarName73 = cell2mat(raw(:, 73));
MazzioJon3.VarName74 = cell2mat(raw(:, 74));
MazzioJon3.VarName75 = cell2mat(raw(:, 75));
MazzioJon3.VarName76 = cell2mat(raw(:, 76));
MazzioJon3.VarName77 = cell2mat(raw(:, 77));
MazzioJon3.VarName78 = cell2mat(raw(:, 78));
MazzioJon3.VarName79 = cell2mat(raw(:, 79));
MazzioJon3.VarName80 = cell2mat(raw(:, 80));
MazzioJon3.VarName81 = cell2mat(raw(:, 81));
MazzioJon3.VarName82 = cell2mat(raw(:, 82));
MazzioJon3.VarName83 = cell2mat(raw(:, 83));
MazzioJon3.VarName84 = cell2mat(raw(:, 84));
MazzioJon3.VarName85 = cell2mat(raw(:, 85));
MazzioJon3.VarName86 = cell2mat(raw(:, 86));
MazzioJon3.VarName87 = cell2mat(raw(:, 87));
MazzioJon3.VarName88 = cell2mat(raw(:, 88));
MazzioJon3.VarName89 = cell2mat(raw(:, 89));
MazzioJon3.VarName90 = cell2mat(raw(:, 90));
MazzioJon3.VarName91 = cell2mat(raw(:, 91));
MazzioJon3.VarName92 = cell2mat(raw(:, 92));
MazzioJon3.VarName93 = cell2mat(raw(:, 93));
MazzioJon3.VarName94 = cell2mat(raw(:, 94));
MazzioJon3.VarName95 = cell2mat(raw(:, 95));
MazzioJon3.VarName96 = cell2mat(raw(:, 96));
MazzioJon3.VarName97 = cell2mat(raw(:, 97));
MazzioJon3.VarName98 = cell2mat(raw(:, 98));
MazzioJon3.VarName99 = cell2mat(raw(:, 99));
MazzioJon3.VarName100 = cell2mat(raw(:, 100));
MazzioJon3.VarName101 = cell2mat(raw(:, 101));
MazzioJon3.VarName102 = cell2mat(raw(:, 102));
MazzioJon3.VarName103 = cell2mat(raw(:, 103));
MazzioJon3.VarName104 = cell2mat(raw(:, 104));
MazzioJon3.VarName105 = cell2mat(raw(:, 105));
MazzioJon3.VarName106 = cell2mat(raw(:, 106));
MazzioJon3.VarName107 = cell2mat(raw(:, 107));
MazzioJon3.VarName108 = cell2mat(raw(:, 108));
MazzioJon3.VarName109 = cell2mat(raw(:, 109));
MazzioJon3.VarName110 = cell2mat(raw(:, 110));
MazzioJon3.VarName111 = cell2mat(raw(:, 111));
MazzioJon3.VarName112 = cell2mat(raw(:, 112));
MazzioJon3.VarName113 = cell2mat(raw(:, 113));
MazzioJon3.VarName114 = cell2mat(raw(:, 114));
MazzioJon3.VarName115 = cell2mat(raw(:, 115));
MazzioJon3.VarName116 = cell2mat(raw(:, 116));
MazzioJon3.VarName117 = cell2mat(raw(:, 117));
MazzioJon3.VarName118 = cell2mat(raw(:, 118));
MazzioJon3.VarName119 = cell2mat(raw(:, 119));
MazzioJon3.VarName120 = cell2mat(raw(:, 120));
MazzioJon3.VarName121 = cell2mat(raw(:, 121));
MazzioJon3.VarName122 = cell2mat(raw(:, 122));
MazzioJon3.VarName123 = cell2mat(raw(:, 123));
MazzioJon3.VarName124 = cell2mat(raw(:, 124));
MazzioJon3.VarName125 = cell2mat(raw(:, 125));
MazzioJon3.VarName126 = cell2mat(raw(:, 126));
MazzioJon3.VarName127 = cell2mat(raw(:, 127));
MazzioJon3.VarName128 = cell2mat(raw(:, 128));
MazzioJon3.VarName129 = cell2mat(raw(:, 129));
MazzioJon3.VarName130 = cell2mat(raw(:, 130));
MazzioJon3.VarName131 = cell2mat(raw(:, 131));
MazzioJon3.VarName132 = cell2mat(raw(:, 132));
MazzioJon3.VarName133 = cell2mat(raw(:, 133));
MazzioJon3.VarName134 = cell2mat(raw(:, 134));
MazzioJon3.VarName135 = cell2mat(raw(:, 135));
MazzioJon3.VarName136 = cell2mat(raw(:, 136));
MazzioJon3.VarName137 = cell2mat(raw(:, 137));
MazzioJon3.VarName138 = cell2mat(raw(:, 138));
MazzioJon3.VarName139 = cell2mat(raw(:, 139));
MazzioJon3.VarName140 = cell2mat(raw(:, 140));
MazzioJon3.VarName141 = cell2mat(raw(:, 141));
MazzioJon3.VarName142 = cell2mat(raw(:, 142));
MazzioJon3.VarName143 = cell2mat(raw(:, 143));
MazzioJon3.VarName144 = cell2mat(raw(:, 144));
MazzioJon3.VarName145 = cell2mat(raw(:, 145));
MazzioJon3.VarName146 = cell2mat(raw(:, 146));
MazzioJon3.VarName147 = cell2mat(raw(:, 147));
MazzioJon3.VarName148 = cell2mat(raw(:, 148));
MazzioJon3.VarName149 = cell2mat(raw(:, 149));
MazzioJon3.VarName150 = cell2mat(raw(:, 150));
MazzioJon3.VarName151 = cell2mat(raw(:, 151));
MazzioJon3.VarName152 = cell2mat(raw(:, 152));
MazzioJon3.VarName153 = cell2mat(raw(:, 153));
MazzioJon3.VarName154 = cell2mat(raw(:, 154));
MazzioJon3.VarName155 = cell2mat(raw(:, 155));
MazzioJon3.VarName156 = cell2mat(raw(:, 156));
MazzioJon3.VarName157 = cell2mat(raw(:, 157));
MazzioJon3.VarName158 = cell2mat(raw(:, 158));
MazzioJon3.VarName159 = cell2mat(raw(:, 159));
MazzioJon3.VarName160 = cell2mat(raw(:, 160));
MazzioJon3.VarName161 = cell2mat(raw(:, 161));
MazzioJon3.VarName162 = cell2mat(raw(:, 162));
MazzioJon3.VarName163 = cell2mat(raw(:, 163));
MazzioJon3.VarName164 = cell2mat(raw(:, 164));
MazzioJon3.VarName165 = cell2mat(raw(:, 165));
MazzioJon3.VarName166 = cell2mat(raw(:, 166));
MazzioJon3.VarName167 = cell2mat(raw(:, 167));
MazzioJon3.VarName168 = cell2mat(raw(:, 168));
MazzioJon3.VarName169 = cell2mat(raw(:, 169));
MazzioJon3.VarName170 = cell2mat(raw(:, 170));
MazzioJon3.VarName171 = cell2mat(raw(:, 171));
MazzioJon3.VarName172 = cell2mat(raw(:, 172));
MazzioJon3.VarName173 = cell2mat(raw(:, 173));
MazzioJon3.VarName174 = cell2mat(raw(:, 174));
MazzioJon3.VarName175 = cell2mat(raw(:, 175));
MazzioJon3.VarName176 = cell2mat(raw(:, 176));
MazzioJon3.VarName177 = cell2mat(raw(:, 177));
MazzioJon3.VarName178 = cell2mat(raw(:, 178));
MazzioJon3.VarName179 = cell2mat(raw(:, 179));
MazzioJon3.VarName180 = cell2mat(raw(:, 180));
MazzioJon3.VarName181 = cell2mat(raw(:, 181));
MazzioJon3.VarName182 = cell2mat(raw(:, 182));
MazzioJon3.VarName183 = cell2mat(raw(:, 183));
MazzioJon3.VarName184 = cell2mat(raw(:, 184));
MazzioJon3.VarName185 = cell2mat(raw(:, 185));
MazzioJon3.VarName186 = cell2mat(raw(:, 186));
MazzioJon3.VarName187 = cell2mat(raw(:, 187));
MazzioJon3.VarName188 = cell2mat(raw(:, 188));
MazzioJon3.VarName189 = cell2mat(raw(:, 189));
MazzioJon3.VarName190 = cell2mat(raw(:, 190));
MazzioJon3.VarName191 = cell2mat(raw(:, 191));
MazzioJon3.VarName192 = cell2mat(raw(:, 192));
MazzioJon3.VarName193 = cell2mat(raw(:, 193));
MazzioJon3.VarName194 = cell2mat(raw(:, 194));
MazzioJon3.VarName195 = cell2mat(raw(:, 195));
MazzioJon3.VarName196 = cell2mat(raw(:, 196));
MazzioJon3.VarName197 = cell2mat(raw(:, 197));
MazzioJon3.VarName198 = cell2mat(raw(:, 198));
MazzioJon3.VarName199 = cell2mat(raw(:, 199));

