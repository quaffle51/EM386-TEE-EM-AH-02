# Octave only supports @CLASS, not classdef
# So use a struct plus functions as a substitute for a utilities class

global octavetex = struct();
octavetex.docdir = pwd();
try
    cd '/home/gordon/Dropbox/OU/M836/TMAs/TMA02/.';
catch
    arg_list = argv()
    if size(arg_list, 1) == 1 && arg_list{1} == '--manual'
    else
        error("Could not find directory .");
    end
end
if dir_in_loadpath(octavetex.docdir)
else
    addpath(octavetex.docdir);
end



octavetex.dependencies = {};
octavetex.created = {};
octavetex._context_raw = '';

function octavetex_formatter(argin)
    disp(argin);
end
octavetex.formatter = @(argin) octavetex_formatter(argin);

function octavetex_before()
end
octavetex.before = @() octavetex_before();

function octavetex_after()
end
octavetex.after = @() octavetex_after();

function octavetex_add_dependencies(varargin)
    global octavetex;
    for i = 1:length(varargin)
        octavetex.dependencies{end+1} = varargin{i};
    end
end
octavetex.add_dependencies = @(varargin) octavetex_add_dependencies(varargin{:});

function octavetex_add_created(varargin)
    global octavetex;
    for i = 1:length(varargin)
        octavetex.created{end+1} = varargin{i};
    end
end
octavetex.add_created = @(varargin) octavetex_add_created(varargin{:});

function octavetex_set_context(argin)
    global octavetex;
    if ~strcmp(argin, octavetex._context_raw)
        octavetex._context_raw = argin;
        hash = struct;
        argin_kv = strsplit(argin, ',');
        for i = 1:length(argin_kv)
            kv = strsplit(argin_kv{i}, '=');
            k = strtrim(kv{1});
            v = strtrim(kv{2});
            hash = setfield(hash, k, v);
        end
        octavetex.context = hash;
    end
end
octavetex.set_context = @(argin) octavetex_set_context(argin);

function out = octavetex_pt_to_in(argin)
    if ischar(argin)
        if length(argin) > 2 && argin(end-1:end) == 'pt'
            out = str2num(argin(1:end-2))/72.27;
        else
            out = str2num(argin)/72.27;
        end
    else
        out = argin/72.27;
    end
end
octavetex.pt_to_in = @(argin) octavetex_pt_to_in(argin);

function out = octavetex_pt_to_cm(argin)
    out = octavetex_pt_to_in(argin)*2.54;
end
octavetex.pt_to_cm = @(argin) octavetex_pt_to_cm(argin);

function out = octavetex_pt_to_mm(argin)
    out = octavetex_pt_to_in(argin)*25.4;
end
octavetex.pt_to_mm = @(argin) octavetex_pt_to_mm(argin);

function out = octavetex_pt_to_bp(argin)
    out = octavetex_pt_to_in(argin)*72;
end
octavetex.pt_to_bp = @(argin) octavetex_pt_to_bp(argin);

function octavetex_cleanup()
    global octavetex;
    fprintf(strcat('=>PYTHONTEX:DEPENDENCIES#', "\n"));
    for i = 1:length(octavetex.dependencies)
        fprintf(strcat(octavetex.dependencies{i}, "\n"));
    end
    fprintf(strcat('=>PYTHONTEX:CREATED#', "\n"));
    for i = 1:length(octavetex.created)
        fprintf(strcat(octavetex.created{i}, "\n"));
    end
end
octavetex.cleanup = @() octavetex_cleanup();

octavetex.id = 'octave_default_0';
octavetex.family = 'octave';
octavetex.session = 'default';
octavetex.restart = '0';

octavetex.command = 'code';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '0';
octavetex.line = '24';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#0#code#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#0#code#', "\n"));
function H = to_standard_form
  n = 8;
  k = 2;
  G=[0, 1, 2, 4, 6, 4, 3, 5; 3, 2, 2, 6, 1, 2, 2, 0];
  G = [mod(3*G(1,:),7); mod(5*G(2,:),7)];
  G = [G(1,:); mod(G(2,:)-G(1,:),7)];
  G = [mod(5*G(1,:),7);mod(1*G(2,:),7)];
  G = [mod(1*G(2,:),7);mod(1*G(1,:),7)];
  A = G(1:2,3:8);
  H = [mod(-A',7),eye(n-k)];
endfunction;

octavetex.after()
octavetex.command = 'code';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '1';
octavetex.line = '38';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#1#code#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#1#code#', "\n"));
function result = syndrome(y, H)
	result = mod(y * H', 7);
endfunction

H = to_standard_form;

y1 = [1, 0, 0, 0, 0, 0, 0, 0];
y2 = [0, 1, 0, 0, 0, 0, 0, 0];
y3 = [0, 0, 1, 0, 0, 0, 0, 0];
y4 = [0, 0, 0, 1, 0, 0, 0, 0];
y5 = [0, 0, 0, 0, 1, 0, 0, 0];
y6 = [0, 0, 0, 0, 0, 1, 0, 0];
y7 = [0, 0, 0, 0, 0, 0, 1, 0];
y8 = [0, 0, 0, 0, 0, 0, 0, 1];

s1 = y1*H';
s2 = y2*H';
s3 = y3*H';
s4 = y4*H';
s5 = y5*H';
s6 = y6*H';
s7 = y7*H';
s8 = y8*H';

y = [4, 5, 6, 3, 2, 0, 3, 6];
s = syndrome(y, H);
octavetex.after()
octavetex.command = 'i';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '2';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#2#i#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#2#i#', "\n"));
disp(disp(y1))
octavetex.after()
octavetex.command = 'i';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '3';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#3#i#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#3#i#', "\n"));
disp(disp(y2))
octavetex.after()
octavetex.command = 'i';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '4';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#4#i#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#4#i#', "\n"));
disp(disp(y8))
octavetex.after()
octavetex.command = 'i';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '5';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#5#i#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#5#i#', "\n"));
disp(disp(y))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '6';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#6#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#6#c#', "\n"));
disp(s)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '7';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#7#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#7#c#', "\n"));
disp(s2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '8';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#8#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#8#c#', "\n"));
disp(mod(s2*3,7))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '9';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#9#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#9#c#', "\n"));
disp(y)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '10';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#10#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#10#c#', "\n"));
disp(y2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '11';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#11#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#11#c#', "\n"));
disp(y2*3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '12';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#12#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#12#c#', "\n"));
disp(y)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '13';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#13#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#13#c#', "\n"));
disp(y2*3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '14';
octavetex.line = '65';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#14#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#14#c#', "\n"));
disp(y-y2*3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '15';
octavetex.line = '69';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#15#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#15#c#', "\n"));
disp(y1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '16';
octavetex.line = '69';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#16#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#16#c#', "\n"));
disp(s1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '17';
octavetex.line = '70';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#17#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#17#c#', "\n"));
disp(y2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '18';
octavetex.line = '70';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#18#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#18#c#', "\n"));
disp(s2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '19';
octavetex.line = '71';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#19#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#19#c#', "\n"));
disp(y3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '20';
octavetex.line = '71';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#20#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#20#c#', "\n"));
disp(s3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '21';
octavetex.line = '72';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#21#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#21#c#', "\n"));
disp(y4)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '22';
octavetex.line = '72';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#22#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#22#c#', "\n"));
disp(s4)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '23';
octavetex.line = '73';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#23#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#23#c#', "\n"));
disp(y5)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '24';
octavetex.line = '73';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#24#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#24#c#', "\n"));
disp(s5)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '25';
octavetex.line = '74';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#25#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#25#c#', "\n"));
disp(y6)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '26';
octavetex.line = '74';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#26#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#26#c#', "\n"));
disp(s6)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '27';
octavetex.line = '75';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#27#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#27#c#', "\n"));
disp(y7)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '28';
octavetex.line = '75';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#28#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#28#c#', "\n"));
disp(s7)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '29';
octavetex.line = '76';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#29#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#29#c#', "\n"));
disp(y8)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '30';
octavetex.line = '76';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#30#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#30#c#', "\n"));
disp(s8)
octavetex.after()
octavetex.command = 'code';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '31';
octavetex.line = '3';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#31#code#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#31#code#', "\n"));
function sy = h_hat(y)

      H =  [0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  0 ;
            0  0  0  1  1  1  1  0  0  0  0  1  1  1  1  0 ;
            0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0 ;
            1  0  1  0  1  0  1  0  1  0  1  0  1  0  1  0 ;
            1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1];

      sy = mod(y * H', 2);
endfunction

y1 = [0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0];
y2 = [0 0 0 1 1 0 0 1 1 1 1 0 0 1 1 1];
y3 = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1];
s1 = h_hat(y1);
s2 = h_hat(y2);
s3 = h_hat(y3);
b1 = mod(s1 * [2^4 2^3 2^2 2^1 2^0]', 2);
b2 = mod(s2 * [2^4 2^3 2^2 2^1 2^0]', 2);
b3 = mod(s3 * [2^4 2^3 2^2 2^1 2^0]', 2);

octavetex.after()
octavetex.command = 'sub';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '32';
octavetex.line = '29';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#32#sub#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#32#sub#', "\n"));
disp("=>PYTHONTEX:FIELD_DELIM#")
disp(strcat(["$\\bm{y} = [$", num2str(y1), "$]$ so $S(\\bm{y}) = \\bm{y}\\hat{H}^T =[$", num2str(sy1=h_hat(y1)), "]"]))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '33';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#33#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#33#c#', "\n"));
disp(sy1(5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '34';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#34#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#34#c#', "\n"));
disp(sy1(1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '35';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#35#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#35#c#', "\n"));
disp(sy1(2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '36';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#36#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#36#c#', "\n"));
disp(sy1(3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '37';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#37#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#37#c#', "\n"));
disp(sy1(4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '38';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#38#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#38#c#', "\n"));
disp([sy1(1) sy1(2) sy1(3) sy1(4)]*[8 4 2 1]')
octavetex.after()
octavetex.command = 'sub';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '39';
octavetex.line = '34';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#39#sub#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#39#sub#', "\n"));
disp("=>PYTHONTEX:FIELD_DELIM#")
disp(strcat(["$\\bm{y} = [$", num2str(y2), "$]$ so $S(\\bm{y}) = \\bm{y}\\hat{H}^T =[$", num2str(sy2=h_hat(y2)), "]"]))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '40';
octavetex.line = '36';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#40#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#40#c#', "\n"));
disp(sy2(5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '41';
octavetex.line = '36';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#41#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#41#c#', "\n"));
disp(sy2(1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '42';
octavetex.line = '36';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#42#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#42#c#', "\n"));
disp(sy2(2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '43';
octavetex.line = '36';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#43#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#43#c#', "\n"));
disp(sy2(3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '44';
octavetex.line = '36';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#44#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#44#c#', "\n"));
disp(sy2(4))
octavetex.after()
octavetex.command = 'sub';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '45';
octavetex.line = '39';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#45#sub#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#45#sub#', "\n"));
disp("=>PYTHONTEX:FIELD_DELIM#")
disp(strcat(["$\\bm{y} = [$", num2str(y3), "$]$ so $S(\\bm{y}) = \\bm{y}\\hat{H}^T =[$", num2str(sy3=h_hat(y3)), "]"]))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '46';
octavetex.line = '41';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#46#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#46#c#', "\n"));
disp(sy3(5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '47';
octavetex.line = '41';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#47#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#47#c#', "\n"));
disp(sy3(1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '48';
octavetex.line = '41';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#48#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#48#c#', "\n"));
disp(sy3(2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '49';
octavetex.line = '41';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#49#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#49#c#', "\n"));
disp(sy3(3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '50';
octavetex.line = '41';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#50#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#50#c#', "\n"));
disp(sy3(4))
octavetex.after()


octavetex.cleanup()
