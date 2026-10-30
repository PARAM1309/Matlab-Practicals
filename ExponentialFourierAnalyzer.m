
function ExponentialFourierAnalyzer

clc;
close all;

A = 1;
T = 10;
w0 = 2*pi/T;
t = linspace(-T/2,T/2,4001);

% GUI WINDOW
fig = uifigure('Name','Exponential Fourier Analyzer',...
    'Position',[50 50 1150 750]);

uilabel(fig,'Text','EXPONENTIAL FOURIER SERIES ANALYZER',...
    'Position',[300 700 600 30],...
    'FontSize',20,'FontWeight','bold');

% INPUT CONTROLS
uilabel(fig,'Text','Waveform:',...
    'Position',[30 655 100 25]);

wave = uidropdown(fig,...
    'Items',{'Even Triangular Wave',...
    'Odd Sawtooth Wave',...
    'Half-Wave Symmetric Square Wave'},...
    'Position',[130 655 250 30]);

uilabel(fig,'Text','Harmonics:',...
    'Position',[430 655 100 25]);

Nfield = uieditfield(fig,'numeric',...
    'Position',[530 655 80 30],...
    'Value',50,'Limits',[1 100],...
    'RoundFractionalValues','on');

uibutton(fig,'Text','ANALYZE',...
    'Position',[700 650 130 40],...
    'ButtonPushedFcn',@analyze);

% PLOT AREAS
ax1 = uiaxes(fig,'Position',[30 390 510 230]);
title(ax1,'Original Waveform');

ax2 = uiaxes(fig,'Position',[600 390 510 230]);
title(ax2,'Fourier Reconstruction');

% COEFFICIENT TABLES
uilabel(fig,'Text','Exponential Fourier Coefficients',...
    'Position',[30 355 400 25],...
    'FontWeight','bold');

tbl1 = uitable(fig,...
    'Position',[30 120 510 230],...
    'ColumnName',{'n','Real(Cn)','Imag(Cn)','Analytical Cn'});

uilabel(fig,'Text','Trigonometric Fourier Coefficients',...
    'Position',[600 355 400 25],...
    'FontWeight','bold');

tbl2 = uitable(fig,...
    'Position',[600 120 510 230],...
    'ColumnName',{'n','an','bn'});

% RESULTS
info = uitextarea(fig,...
    'Position',[30 15 1080 80],...
    'Editable','off');

analyze();

% =====================================
% ANALYSIS FUNCTION
% =====================================

    function analyze(~,~)

        N = round(Nfield.Value);
        choice = wave.Value;

        % SIGNAL GENERATION
        switch choice

            case 'Even Triangular Wave'
                x = A*(1-4*abs(t)/T);
                symmetry = 'Even';
                shortcut = 'bn=0, even harmonics=0';
                reduction = 75;

            case 'Odd Sawtooth Wave'
                x = 2*A*t/T;
                symmetry = 'Odd';
                shortcut = 'an=0';
                reduction = 50;

            otherwise
                x = A*sign(sin(w0*t));
                symmetry = 'Odd and Half-wave';
                shortcut = 'an=0, even harmonics=0';
                reduction = 75;

        end

        % EXPONENTIAL FOURIER COEFFICIENTS
        n = (-N:N)';
        Cn = zeros(size(n));
        Can = NaN(size(n));

        for k = 1:length(n)
            Cn(k) = trapz(t,...
                x.*exp(-1j*n(k)*w0*t))/T;
        end

        % ANALYTICAL TRIANGULAR COEFFICIENTS
        if strcmp(choice,'Even Triangular Wave')

            for k = 1:length(n)
                if n(k) ~= 0 && mod(abs(n(k)),2)==1
                    Can(k) = 4*A/(pi^2*n(k)^2);
                else
                    Can(k) = 0;
                end
            end

            err = max(abs(Cn-Can));
            errorText = sprintf('Maximum error = %.8f',err);

        else
            errorText = 'Analytical comparison: triangular wave only';
        end

        % TRIGONOMETRIC FOURIER COEFFICIENTS
        a0 = 2*trapz(t,x)/T;
        an = zeros(N,1);
        bn = zeros(N,1);

        for k = 1:N

            if strcmp(choice,'Even Triangular Wave')

                if mod(k,2)==1
                    an(k) = 2*trapz(t,...
                        x.*cos(k*w0*t))/T;
                end

            elseif strcmp(choice,'Odd Sawtooth Wave')

                bn(k) = 2*trapz(t,...
                    x.*sin(k*w0*t))/T;

            else

                if mod(k,2)==1
                    bn(k) = 2*trapz(t,...
                        x.*sin(k*w0*t))/T;
                end

            end
        end

        % FOURIER RECONSTRUCTION
        xr = (a0/2)*ones(size(t));

        for k = 1:N
            xr = xr + an(k)*cos(k*w0*t)...
                + bn(k)*sin(k*w0*t);
        end

        % ORIGINAL SIGNAL
        plot(ax1,t,x,'b','LineWidth',1.8);
        grid(ax1,'on');
        xlabel(ax1,'Time (seconds)');
        ylabel(ax1,'Amplitude');
        title(ax1,choice);

        % RECONSTRUCTED SIGNAL
        plot(ax2,t,x,'k--','LineWidth',1.2);
        hold(ax2,'on');
        plot(ax2,t,xr,'r','LineWidth',1.8);
        hold(ax2,'off');

        grid(ax2,'on');
        xlabel(ax2,'Time (seconds)');
        ylabel(ax2,'Amplitude');
        legend(ax2,{'Original','Reconstructed'});
        title(ax2,'Fourier Reconstruction');

        % DISPLAY TABLES
        tbl1.Data = [n,real(Cn),imag(Cn),Can];
        tbl2.Data = [(1:N)',an,bn];

        % SYMMETRY AND RESULTS
        info.Value = {
            ['Time Period = ',num2str(T),' seconds']
            ['Amplitude = ',num2str(A),...
            ' | Harmonics = ',num2str(N),...
            ' | Frequency = ',num2str(w0),' rad/sec']
            ['Symmetry = ',symmetry,...
            ' | Shortcut = ',shortcut]
            ['Calculation reduction = ',...
            num2str(reduction),'%']
            errorText
        };

    end

end
