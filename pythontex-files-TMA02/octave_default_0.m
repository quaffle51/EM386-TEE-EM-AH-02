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
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '51';
octavetex.line = '25';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#51#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#51#c#', "\n"));
disp((5-1)/2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '52';
octavetex.line = '25';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#52#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#52#c#', "\n"));
disp((5-1)/2)
octavetex.after()
octavetex.command = 'code';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '53';
octavetex.line = '3';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#53#code#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#53#code#', "\n"));
	function i = Inverse(x)
	  T = [
	      0 0 0 0 0 0 0;
	      0 1 2 3 4 5 6;
	      0 2 4 6 1 3 5;
	      0 3 6 2 5 1 4;
	      0 4 1 5 2 6 3;
	      0 5 3 1 6 4 2;
	      0 6 5 4 3 2 1
	    ];
	  i = 0;	
	  if x != 0
	    for j = [1 2 3 4 5 6 7]
	      if T(x+1,j) != 1
	        i++;
	      else
	        break;
	      end
	    end
	  end
	
	endfunction
	
	y = [3 2 4 6 6 4];
	
	
	H = [1^0 2^0 3^0 4^0 5^0 6^0;
	     1^1 2^1 3^1 4^1 5^1 6^1;
	     1^2 2^2 3^2 4^2 5^2 6^2;
	     1^3 2^3 3^3 4^3 5^3 6^3
	    ];
	
	H  = mod(H, 7);
	
	y1 = y(1);
	y2 = y(2);
	y3 = y(3);
	y4 = y(4);
	y5 = y(5);
	y6 = y(6);
	
	S = mod(y*H', 7);
	
	S1a = S(1);
	S2a = S(2);
	S3a = S(3);
	S4a = S(4);
	H11 = H(1,1);
	
	function result =  Sj(j, y)
		result = 0;
		for i=[1 2 3 4 5 6]
			result += y(i) * i^(j-1);
		end
		result = mod(result, 7);
	endfunction
	
	S1 = Sj(1,y);
	S2 = Sj(2,y);
	S3 = Sj(3,y);
	S4 = Sj(4,y);
	
	A1 = 4;
	A2 = 1;
	B1 = 4;
	B2 = 2;
	
	#factor (1 + 4x + 2x^2) mod 7
	Z2 = mod(-3,7);
	Z1 = mod(-6,7);
	z1 = Inverse(Z1);
	z2 = Inverse(Z2);
	
	m1 = 2;
	m2 = 2;
	X1 = 1;
	X2 = 2;
	
	v = [m1 m2 4 6 6 4];
	sy = mod(v * H', 7);
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '54';
octavetex.line = '84';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#54#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#54#c#', "\n"));
disp(y)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '55';
octavetex.line = '87';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#55#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#55#c#', "\n"));
disp(y)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '56';
octavetex.line = '112';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#56#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#56#c#', "\n"));
disp(S1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '57';
octavetex.line = '112';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#57#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#57#c#', "\n"));
disp(S2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '58';
octavetex.line = '112';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#58#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#58#c#', "\n"));
disp(S3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '59';
octavetex.line = '112';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#59#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#59#c#', "\n"));
disp(S4)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '60';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#60#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#60#c#', "\n"));
disp(S1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '61';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#61#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#61#c#', "\n"));
disp(S2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '62';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#62#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#62#c#', "\n"));
disp(S1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '63';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#63#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#63#c#', "\n"));
disp(S3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '64';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#64#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#64#c#', "\n"));
disp(S2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '65';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#65#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#65#c#', "\n"));
disp(S1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '66';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#66#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#66#c#', "\n"));
disp(S4)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '67';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#67#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#67#c#', "\n"));
disp(S3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '68';
octavetex.line = '130';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#68#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#68#c#', "\n"));
disp(S2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '69';
octavetex.line = '154';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#69#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#69#c#', "\n"));
disp(B1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '70';
octavetex.line = '154';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#70#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#70#c#', "\n"));
disp(4 * B1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '71';
octavetex.line = '154';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#71#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#71#c#', "\n"));
disp(6 + 16)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '72';
octavetex.line = '154';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#72#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#72#c#', "\n"));
disp(mod(6 + 16,7))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '73';
octavetex.line = '155';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#73#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#73#c#', "\n"));
disp(A1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '74';
octavetex.line = '155';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#74#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#74#c#', "\n"));
disp(A2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '75';
octavetex.line = '155';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#75#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#75#c#', "\n"));
disp(B1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '76';
octavetex.line = '155';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#76#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#76#c#', "\n"));
disp(B2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '77';
octavetex.line = '160';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#77#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#77#c#', "\n"));
disp(A1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '78';
octavetex.line = '160';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#78#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#78#c#', "\n"));
disp(B1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '79';
octavetex.line = '160';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#79#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#79#c#', "\n"));
disp(B2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '80';
octavetex.line = '161';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#80#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#80#c#', "\n"));
disp(Z1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '81';
octavetex.line = '161';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#81#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#81#c#', "\n"));
disp(Z2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '82';
octavetex.line = '161';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#82#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#82#c#', "\n"));
disp(z1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '83';
octavetex.line = '161';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#83#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#83#c#', "\n"));
disp(z2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '84';
octavetex.line = '170';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#84#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#84#c#', "\n"));
disp(z2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '85';
octavetex.line = '170';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#85#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#85#c#', "\n"));
disp(z2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '86';
octavetex.line = '179';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#86#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#86#c#', "\n"));
disp(z2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '87';
octavetex.line = '179';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#87#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#87#c#', "\n"));
disp(z2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '88';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#88#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#88#c#', "\n"));
disp(mod(1 + 1 + 0 + 5 + 2 + 5,7))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '89';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#89#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#89#c#', "\n"));
disp(mod( 1 + 2 + 3*0 + 4*5 + 5*2 + 6*5,7))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '90';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#90#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#90#c#', "\n"));
disp(mod(1 + 4 + 9*0 + 16*5 + 25*2 + 36*5,7))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '91';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#91#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#91#c#', "\n"));
disp(mod(1 + 8 + 27*0 + 64*5 + 125*2 + 216*5,7))
octavetex.after()
octavetex.command = 'code';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '92';
octavetex.line = '3';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#92#code#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#92#code#', "\n"));
	q = 7;
	y1 = [1 1 0 5 2 5];
	y2 = [1 0 4 6 6 4];
	A = [4 6 6 4;3 6 3 1];
	a11 = A(1,1);
	a12 = A(1,2);
	a13 = A(1,3);
	a14 = A(1,4);
	
	a21 = A(2,1);
	a22 = A(2,2);
	a23 = A(2,3);
	a24 = A(2,4);
	G = [eye(2),A];
	
	v = [4, 5]
	test1 = mod(v*G, q);
	GT = G';
	check = mod(H*G', q);
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '93';
octavetex.line = '31';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#93#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#93#c#', "\n"));
disp(y)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '94';
octavetex.line = '46';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#94#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#94#c#', "\n"));
disp(y(1,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '95';
octavetex.line = '46';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#95#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#95#c#', "\n"));
disp(y(1,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '96';
octavetex.line = '46';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#96#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#96#c#', "\n"));
disp(y(1,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '97';
octavetex.line = '46';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#97#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#97#c#', "\n"));
disp(y(1,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '98';
octavetex.line = '56';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#98#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#98#c#', "\n"));
disp(a11)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '99';
octavetex.line = '57';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#99#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#99#c#', "\n"));
disp(a12)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '100';
octavetex.line = '58';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#100#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#100#c#', "\n"));
disp(a13)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '101';
octavetex.line = '59';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#101#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#101#c#', "\n"));
disp(a14)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '102';
octavetex.line = '60';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#102#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#102#c#', "\n"));
disp(a21)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '103';
octavetex.line = '61';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#103#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#103#c#', "\n"));
disp(a22)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '104';
octavetex.line = '62';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#104#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#104#c#', "\n"));
disp(a23)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '105';
octavetex.line = '63';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#105#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#105#c#', "\n"));
disp(a24)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '106';
octavetex.line = '67';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#106#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#106#c#', "\n"));
disp(G(1,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '107';
octavetex.line = '68';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#107#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#107#c#', "\n"));
disp(G(1,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '108';
octavetex.line = '69';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#108#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#108#c#', "\n"));
disp(G(1,3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '109';
octavetex.line = '70';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#109#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#109#c#', "\n"));
disp(G(1,4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '110';
octavetex.line = '71';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#110#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#110#c#', "\n"));
disp(G(1,5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '111';
octavetex.line = '72';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#111#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#111#c#', "\n"));
disp(G(1,6))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '112';
octavetex.line = '73';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#112#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#112#c#', "\n"));
disp(G(2,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '113';
octavetex.line = '74';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#113#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#113#c#', "\n"));
disp(G(2,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '114';
octavetex.line = '75';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#114#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#114#c#', "\n"));
disp(G(2,3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '115';
octavetex.line = '76';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#115#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#115#c#', "\n"));
disp(G(2,4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '116';
octavetex.line = '77';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#116#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#116#c#', "\n"));
disp(G(2,5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '117';
octavetex.line = '78';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#117#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#117#c#', "\n"));
disp(G(2,6))
octavetex.after()
octavetex.command = 'code';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '118';
octavetex.line = '3';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#118#code#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#118#code#', "\n"));
	v1 = [3 2 4 1 3 0];
	c1 = mod([v1(1,1) v1(1,2)]*G, q);
	v2 = [4 5 2 0 6 6];
	c2 = mod([v2(1,1) v2(1,2)]*G, q);
	c3 = mod([2 5]*G, q);
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '119';
octavetex.line = '9';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#119#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#119#c#', "\n"));
disp(v1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '120';
octavetex.line = '11';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#120#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#120#c#', "\n"));
disp(v1(1,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '121';
octavetex.line = '11';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#121#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#121#c#', "\n"));
disp(v1(1,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '122';
octavetex.line = '13';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#122#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#122#c#', "\n"));
disp(G(1,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '123';
octavetex.line = '14';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#123#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#123#c#', "\n"));
disp(G(1,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '124';
octavetex.line = '15';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#124#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#124#c#', "\n"));
disp(G(1,3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '125';
octavetex.line = '16';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#125#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#125#c#', "\n"));
disp(G(1,4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '126';
octavetex.line = '17';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#126#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#126#c#', "\n"));
disp(G(1,5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '127';
octavetex.line = '18';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#127#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#127#c#', "\n"));
disp(G(1,6))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '128';
octavetex.line = '19';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#128#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#128#c#', "\n"));
disp(G(2,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '129';
octavetex.line = '20';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#129#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#129#c#', "\n"));
disp(G(2,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '130';
octavetex.line = '21';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#130#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#130#c#', "\n"));
disp(G(2,3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '131';
octavetex.line = '22';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#131#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#131#c#', "\n"));
disp(G(2,4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '132';
octavetex.line = '23';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#132#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#132#c#', "\n"));
disp(G(2,5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '133';
octavetex.line = '24';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#133#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#133#c#', "\n"));
disp(G(2,6))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '134';
octavetex.line = '26';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#134#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#134#c#', "\n"));
disp(c1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '135';
octavetex.line = '28';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#135#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#135#c#', "\n"));
disp(v1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '136';
octavetex.line = '28';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#136#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#136#c#', "\n"));
disp(c1)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '137';
octavetex.line = '30';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#137#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#137#c#', "\n"));
disp(v2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '138';
octavetex.line = '32';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#138#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#138#c#', "\n"));
disp(v2(1,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '139';
octavetex.line = '32';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#139#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#139#c#', "\n"));
disp(v2(1,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '140';
octavetex.line = '34';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#140#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#140#c#', "\n"));
disp(G(1,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '141';
octavetex.line = '35';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#141#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#141#c#', "\n"));
disp(G(1,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '142';
octavetex.line = '36';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#142#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#142#c#', "\n"));
disp(G(1,3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '143';
octavetex.line = '37';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#143#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#143#c#', "\n"));
disp(G(1,4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '144';
octavetex.line = '38';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#144#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#144#c#', "\n"));
disp(G(1,5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '145';
octavetex.line = '39';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#145#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#145#c#', "\n"));
disp(G(1,6))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '146';
octavetex.line = '40';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#146#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#146#c#', "\n"));
disp(G(2,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '147';
octavetex.line = '41';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#147#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#147#c#', "\n"));
disp(G(2,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '148';
octavetex.line = '42';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#148#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#148#c#', "\n"));
disp(G(2,3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '149';
octavetex.line = '43';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#149#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#149#c#', "\n"));
disp(G(2,4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '150';
octavetex.line = '44';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#150#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#150#c#', "\n"));
disp(G(2,5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '151';
octavetex.line = '45';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#151#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#151#c#', "\n"));
disp(G(2,6))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '152';
octavetex.line = '47';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#152#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#152#c#', "\n"));
disp(c2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '153';
octavetex.line = '49';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#153#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#153#c#', "\n"));
disp(v2)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '154';
octavetex.line = '53';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#154#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#154#c#', "\n"));
disp(G(1,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '155';
octavetex.line = '54';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#155#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#155#c#', "\n"));
disp(G(1,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '156';
octavetex.line = '55';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#156#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#156#c#', "\n"));
disp(G(1,3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '157';
octavetex.line = '56';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#157#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#157#c#', "\n"));
disp(G(1,4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '158';
octavetex.line = '57';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#158#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#158#c#', "\n"));
disp(G(1,5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '159';
octavetex.line = '58';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#159#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#159#c#', "\n"));
disp(G(1,6))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '160';
octavetex.line = '59';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#160#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#160#c#', "\n"));
disp(G(2,1))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '161';
octavetex.line = '60';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#161#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#161#c#', "\n"));
disp(G(2,2))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '162';
octavetex.line = '61';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#162#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#162#c#', "\n"));
disp(G(2,3))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '163';
octavetex.line = '62';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#163#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#163#c#', "\n"));
disp(G(2,4))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '164';
octavetex.line = '63';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#164#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#164#c#', "\n"));
disp(G(2,5))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '165';
octavetex.line = '64';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#165#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#165#c#', "\n"));
disp(G(2,6))
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '166';
octavetex.line = '66';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#166#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#166#c#', "\n"));
disp(c3)
octavetex.after()
octavetex.command = 'c';
octavetex.set_context('');
octavetex.args = '';
octavetex.instance = '167';
octavetex.line = '68';

octavetex.before()

fprintf(strcat('=>PYTHONTEX:STDOUT#167#c#', "\n"));
fprintf(stderr, strcat('=>PYTHONTEX:STDERR#167#c#', "\n"));
disp(c3)
octavetex.after()


octavetex.cleanup()
