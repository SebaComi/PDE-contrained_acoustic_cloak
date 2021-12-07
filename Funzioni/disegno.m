function [] = disegno(vertices,elements,nodi_domini,Y,Y_inc,salva,nome_video,titolo)
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex');
h = figure(1);
if ~exist('salva','var')
    salva = false;
elseif salva
%     set(gcf, 'Position', [500 2 660 780])
    set(gcf,'Color','w')
%     nome_video = 'Barca - No cloak';
    frame_per_sec = 12;
%     time_scale = 1;
%     posizione_fig = get(gcf,'Position');
    caricam = waitbar(0,'Let''s create the video!');
end
GIF = true;
if ~exist('titolo','var')
    griglia = false;
else
    griglia = true;
end

hold on
% set(gcf,'Position',[ 100  400  1250 350])
bordi = 'off';


Y_somma = sum([Y,Y_inc],2);
Y_somma(~nodi_domini(:,1) & ~nodi_domini(:,2)) = 1i*NaN;
% massim = max(max(abs([Y, Y_somma, Y_inc])));
massim = 1;
Y(~nodi_domini(:,1) & ~nodi_domini(:,2),:) = 1i*NaN;
Y_inc(~nodi_domini(:,1) & ~nodi_domini(:,2)) = 1i*NaN;
TIME = 20;
fattore = 1.4/.8     *0.5;
% figure(1)
figura_tipica([],15/fattore)
for ii = 1:TIME+1
    aa = 2*pi/TIME * (ii-1);
    clf
% %     subplot(131); hold off
set(gcf,'Position',[489 343 [560 420]/fattore])
    pdeplot(vertices,elements(1:3,:),'xydata',real(Y_somma)*sin(aa) + imag(Y_somma)*cos(aa),'xystyle','interp',...
            'zdata',real(Y)*sin(aa) + imag(Y)*cos(aa),'zstyle','continuous','colorbar','on','mesh',bordi);
%     xlabel('x'); ylabel('y'); title('Scattered wave')
xlim([min(vertices(1,:)), max(vertices(1,:))])
ylim([min(vertices(2,:)), max(vertices(2,:))])

set(gcf,'Position',[489 343 [560 420]/fattore])
    caxis([-massim,massim]); colormap(jet); view(2); box on; axis equal %square;% axis equal
    if ~griglia
        axis off;    box off;
        c = colorbar('FontSize',15/fattore,'TickLabelInterpreter','latex');
        if ii < 2
            x = get(c,'Position');
            x1 = get(gca,'position');
            x(3) = x(3)/fattore;
        end
        set(c,'Position',x);
        set(gca,'position',x1);
    else
        title(titolo)
%         xlim([-1 1]*1.05*0.626); ylim([-1 1]*1.05*0.626)
        xlim([-1 1]*1.05); ylim([-1 1]*1.05)
    figura_tipica([],15/fattore)

        if ii < 2
%             x = get(c,'Position');
%             x(3) = x(3)/fattore;
            x1 = get(gca,'position');
        end
% %         set(c,'Position',x);
        set(gca,'position',x1);
    end
    colorbar off
%     subplot(132);
%     pdeplot(vertices,elements(1:3,:),'xydata',real(Y_inc)*sin(aa) + imag(Y_inc)*cos(aa),'xystyle','interp',...
%             'zdata',real(Y_inc)*sin(aa) + imag(Y_inc)*cos(aa),'zstyle','continuous','colorbar','on','mesh',bordi);
%     xlabel('x'); ylabel('y'); title('Incident wave')
%     caxis([-massim,massim]); colormap(jet); view(2);  box on; axis square
%    
%     
%     subplot(133);
%     pdeplot(vertices,elements(1:3,:),'xydata',real(Y_somma)*sin(aa) + imag(Y_somma)*cos(aa),'xystyle','interp',...
%             'zdata',real(Y_somma)*sin(aa) + imag(Y_somma)*cos(aa),'zstyle','continuous','colorbar','on','mesh',bordi);
%     xlabel('x'); ylabel('y'); title('Total wave')
%     caxis([-massim,massim]); colormap(jet); view(2);  box on; axis square
    
    drawnow
    if salva
    %   Se si fa il video
        if GIF
            % Capture the plot as an image 
            frame = getframe(h); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            % Write to the GIF File 
            if ii == 1 
                imwrite(imind,cm,[nome_video, '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1); 
            else 
                imwrite(imind,cm,[nome_video, '.gif'],'gif','WriteMode','append','DelayTime',0.1); 
            end 
        else
            F(ii) = getframe(gcf);
        end
        waitbar(ii/TIME,caricam,'Acquisition of frames')
    end
    
end


if salva && ~GIF
    waitbar(1,caricam,'Acquisition completed')
%% create the video writer
F(ii+1:3) = F(end);
nome_video = [nome_video, '.mj2']; % .mp4

writerObj = VideoWriter(nome_video,'Motion JPEG 2000'); % 'MPEG-4'
writerObj.FrameRate = frame_per_sec;         % set images per second

% open the video writer
open(writerObj);
% write the frames to the video
TOT = length(F);

for i=1:TOT
    % convert the image to a frame
    frame = F(i);    
    writeVideo(writerObj, frame);
    waitbar(i/TOT,caricam,'Creating video...')
end
% close the writer object
close(writerObj);

waitbar(1,caricam,'Video created')
beep
pause(2)
close(caricam)
close all

disp([nome_video, ':    created'])
end

end
%% Nested

function figura_tipica(polare,font_size)
    fattore = 1.4/.8;
    if ~exist('polare','var'),  polare = false;     end
    if isempty(polare),         polare = false;     end
    if ~exist('font_size','var'), font_size = 13;   end
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(0,'defaultTextInterpreter','latex');
    
    set(gcf,'Color','w')
    set(gca,'FontSize',font_size)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontName','Computer Modern')

    if ~polare
        set(gcf,'Position',[489 343 [560 420]/fattore])
        hhh = get(gca,'XLabel');
        hhh.Interpreter = 'latex';
        hhh = get(gca,'YLabel');
        hhh.Interpreter = 'latex';
    else
        set(gcf,'Position',[489 343 [560 420]/fattore])
    end
    grid on;     box on


    hhh = findobj(gcf, 'Type', 'Legend');
    if ~isempty(hhh)
        set(hhh,'FontName','Computer Modern');
        set(hhh,'Interpreter','latex')
    end
    
    hhh = get(gca, 'title');
    if ~isempty(hhh)
        set(hhh,'FontName','Computer Modern');
        set(hhh,'Interpreter','latex')
    end
%     
%     hhh = get(gca,'Legend');
%     if ~isempty(hhh)
%         set(hhh,'Interpreter','latex')
%     end
    x1=get(gca,'position');
    C = findall(gcf,'type','ColorBar');
    for ii=1:length(C)
        
        C(ii).TickLabelInterpreter = 'latex';
        x=get(C(ii),'Position');
        set(C(ii),'Position',[x(1:2) x(3)/1.4 x(4)])
        set(gca,'position',x1);
    end
end

