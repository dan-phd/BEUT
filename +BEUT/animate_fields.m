function animate_fields(dim, varargin)
%Animate 1D or 2D fields in time with media controls.
%
% animate_fields(dim,'parameter',parameter_value) plots any number of
% 1D/2D fields(x,k) side by side where x is in space and k is in time.
% The present time step shown onscreen can be manipulated using play/stop
% and slider buttons. dim specifies the working spatial dimension (1 or 2).
% 
% Use the following parameters to apply to all fields:
%  'domain' to manually set the domain
%  'delay' to set a delay in seconds between each frame [0]
%  'save' followed by video file location to save the movie
%  'fps' to set the playback frame rate of the saved movie [20]
%  'fullscreen' to play the movie in fullscreen
%  'min_amplitude' to set the minimum y-axis
%  'max_amplitude' to set the maximum y-axis
%  'grid' to set the grid 'on' or 'off' ['off']
%  'dimensions' to specify plotting in 2 or 3 dimensions [3]
%  'colormap' to define colormap ['jet']
%  'skipTimesteps' to skip a number of timesteps for each frame [1]
% 
% Use the following parameters after the field to apply only to that field:
%  'min_amplitude' to set the minimum y-axis of the plot
%  'max_amplitude' to set the maximum y-axis of the plot
%  'grid' to set the grid 'on' or 'off' ['off']
%  'dimensions' to specify plotting in 2 or 3 dimensions [2]
%  'overlay' to set coordinates of vertices to overlay on the field
% 
% Any other 'parameter' will be a field name and should be followed by
% its associated field values.
% 
% Example
%   NT=100;
%   time = linspace(0,20*pi,NT);
%   domain = (0:0.01:4*pi)';
%   for k=1:NT
%       sine(:,k) = sin(domain-time(k));
%       cosine(:,k) = cos(domain-time(k));
%       tangent(:,k) = tan(domain-time(k));
%       secant(:,k) = sec(domain-time(k));
%   end
%   animate_fields(1,'domain',domain, 'grid','on', 'save','file.avi', 'fullscreen',...
%       'Sine',sine,...
%       'Cosine',cosine, 'overlay',[5 0; 10 0],...
%       'Tangent',tangent, 'min_amplitude',-1, 'max_amplitude',1, 'grid','off',...
%       'Secant',secant, 'min_amplitude',-2, 'max_amplitude',2)
% 
%   will animate the 1D trig functions side by side in fullscreen with
%   grids showing for all but the tangent, then save the movie to file.avi.
%
%   Note that all fields must be the same size, where size(sine,2) is the
%   number of time steps and size(sine,1) equals size(domain,1).
%

%   Author: Daniel Simmons - DansPhD.com
%   Edited: 23/10/2015


% Default plot types
switch dim
    case 1, default_plot_type='plot';
    case 2, default_plot_type='surf';           % can also support quiver and quiver3
    otherwise, error('dim must be 1 or 2')      % 3D not supported yet
end

% Find values for parameter inputs
show_grid = 'off'; frame_delay = 0; min_amp_given=0; max_amp_given=0; domain_x=0; domain_y=0;
save_location=0; fps=20; fullscreen=0; dimensions=2; customcolor='jet'; skipTimesteps = 1;
global stop_please k;
nArgs = length(varargin); i=0; a=1;
while a < nArgs+1
    if ischar(varargin{a})
        switch varargin{a}
            % Global parameters
            case 'delay'
                frame_delay=varargin{a+1};
                a=a+2;
            case 'save'
                save_location=varargin{a+1};
                a=a+2;
            case 'fps'
                fps=varargin{a+1};
                a=a+2;
            case 'colormap'
                customcolor=varargin{a+1};
                a=a+2;
            case 'fullscreen'
                fullscreen=1;
                a=a+1;
            case 'skipTimesteps'
                skipTimesteps=varargin{a+1};
                a=a+2;
            case 'domain'
                domain_x=varargin{a+1};
                a=a+1;
                if dim>=2
                    domain_y=varargin{a+1};
                    a=a+1;
                end
                if dim==3
                    domain_z=varargin{a+1};
                    a=a+1;
                end
                
            % Parameters specific to certain field
            case 'PlotType'
                if i==0
                    error('plot_type must be specified after each field');
                else
                    field(i).plot_type=varargin{a+1};
                    if strcmp(varargin{a+1},'patch')
                        field(i).custom_plot.type = 'patch';
                        a=a+2;
                    end
                end
                a=a+2;
            case 'grid'
                if i==0
                    show_grid=varargin{a+1};
                else
                    field(i).show_grid=varargin{a+1};
                end
                a=a+2;
            case 'dimensions'
                if i==0
                    dimensions=varargin{a+1};
                else
                    field(i).dimensions=varargin{a+1};
                end
                a=a+2;
            case 'overlay'
                if i==0
                    error('overlay can only be specified after the relevant field')
                else
                    field(i).overlay=varargin{a+1};
                    field(i).plot_overlay = true;

                    assert(size(field(i).overlay,2)==2,...
                        ['overlay must be a matrix of coordinates '...
                        'e.g. [0 0; 1 0; 1 1; 0 1]'])
                end
                a=a+2;
            case 'max_amplitude'
                if i==0
                    max_amp=varargin{a+1};
                    max_amp_given=1;
                else
                    field(i).max_amp=varargin{a+1};
                    field(i).max_amp_given=1;
                end
                a=a+2;
            case 'min_amplitude'
                if i==0
                    min_amp=varargin{a+1};
                    min_amp_given=1;
                else
                    field(i).min_amp=varargin{a+1};
                    field(i).min_amp_given=1;
                end
                a=a+2;
            
            otherwise
                
                % Create a struct for each field
                i=i+1;
                field(i).name = varargin{a};
                field(i).value = varargin{a+1};
                
                % Assign parameters
                field(i).plot_type = default_plot_type;
                field(i).vector_field = false;
                field(i).show_grid=show_grid;
                field(i).dimensions=dimensions;
                field(i).plot_overlay = false;
                field(i).max_amp_given=max_amp_given;
                if max_amp_given~=0, field(i).max_amp=max_amp; end
                field(i).min_amp_given=min_amp_given;
                if min_amp_given~=0, field(i).min_amp=min_amp; end
                
                % make sure the fields are all the same size
                if ~isstruct(field(i).value)     % if not vector fields (magnitudes only)
                    
                    assert(size(field(i).value,1)==size(field(1).value,1),...
                        'Fields must all be the same size. Field ''%s'' is not the same size as field ''%s''',...
                        field(i).name,field(1).name)
                    
                end
                
                a=a+2;
                
        end
    else
        a=a+1;
    end
end
                    

% Total fields
num_fields = i;

% make sure all vector fields are the same size and plot types are set
for i=1:num_fields
    
    if isstruct(field(i).value)     % if vector fields or patch data
        field_names = fieldnames(field(i).value);
        
        % if no plot type was given, make assumption based on input
        if strcmp(field(i).plot_type,'plot')||strcmp(field(i).plot_type,'surf')
            if any(strcmp(field_names,'u'))
                if any(strcmp(field_names,'w'))
                    field(i).plot_type = 'quiver3';
                else
                    field(i).plot_type = 'quiver';
                end
            else
                % the only other plot option is patch
                field(i).plot_type = 'patch';
            end
        end
        
        
        if strcmp(field(i).plot_type,'quiver')
            assert(numel(field_names)==2,...
                'The 2D vector field ''%s'' must have 2 components (u and v)',field(i).name)
            check_for_fieldname('u',field_names,field(i).name);
            check_for_fieldname('v',field_names,field(i).name);
            
            assert(size(field(i).value.v,1)==size(field(1).value.u,1),...
                'Component ''v'' of field ''%s'' is not the same size as component ''u''',...
                field(i).name)
            
            % So we don't have to check each time if the field value is a struct, we make the
            % field "value" equal to one of its components
            field(i).u = field(i).value.u;
            field(i).v = field(i).value.v;
            field(i).value = field(i).u;
            
            field(i).vector_field=true;
            
        elseif strcmp(field(i).plot_type,'quiver3')
            
            assert(numel(field_names)==3,...
                'The 3D vector field ''%s'' must have 3 components (u, v and w)',field(i).name)
            check_for_fieldname('u',field_names,field(i).name);
            check_for_fieldname('v',field_names,field(i).name);
            check_for_fieldname('w',field_names,field(i).name);
            
            assert(size(field(i).value.v,1)==size(field(1).value.u,1),...
                'Component ''v'' of field ''%s'' is not the same size as component ''u''',...
                field(i).name)
            assert(size(field(i).value.v,1)==size(field(1).value.u,1),...
                'Component ''w'' of field ''%s'' is not the same size as component ''u''',...
                field(i).name)
            
            field(i).u = field(i).value.u;
            field(i).v = field(i).value.v;
            field(i).w = field(i).value.w;
            field(i).value = field(i).u;
            
            field(i).vector_field=true;
            
        elseif strcmp(field(i).plot_type,'patch')
            
            check_for_fieldname('Vertices',field_names,field(i).name);
            check_for_fieldname('Faces',field_names,field(i).name);
            check_for_fieldname('FaceVertexCData',field_names,field(i).name);
            
            assert( size(field(i).value.Faces,1)==size(field(i).value.FaceVertexCData,1) | ...
                    size(field(i).value.Vertices,1)==size(field(i).value.FaceVertexCData,1),...
                ['The data contained in ''FaceVertexCData'' of field ''%s'' does not match the'...
                'correspondng number of ''Faces'' or number of ''Vertices'''],...
                field(i).name)
            
            % save the Z data as value and change the dim so that correct number of
            % timesteps and domain size can be calculated further on
            field(i).patches = field(i).value;
            field(i).value = field(i).value.FaceVertexCData;
            dim = 1;
            
        else
            error('Only quiver, quiver3 and patch support struct inputs')
            
        end
    end
end

    function check_for_fieldname(name,struct_titles,field_title)
        assert(any(strcmp(struct_titles,name)),...
            'Field ''%s'' must contain ''%s'' in struct',field_title,name)
    end
    

% Total time steps
switch dim
    case 1, NT = size(field(1).value,2);
    case 2, NT = size(field(1).value,3);
    case 3, NT = size(field(1).value,4);
end


% Check domain, or set default domain if one not given
switch dim
    case 1
        
        if any(domain_x~=0)
            assert(numel(domain_x)==size(field(1).value,1), ...
                'Domain must be the same size as number of field columns')
        else
            domain_x = 1:size(field(1).value,1);
        end
        
    case 2
        
        if any(any(domain_x~=0))
            assert(size(domain_x,1)==size(field(1).value,1), ...
                'Domain must be the same size as field (in a single timestep)')
        else
            [domain_x,domain_y] = deal(1:size(field(1).value,1),1:size(field(1).value,2));
        end
        
    case 3
        
        if any(any(any(domain_x~=0)))
            assert(size(domain_x,1)==size(field(1).value,1), ...
                'Domain must be the same size as field (in a single timestep)')
        else
            [domain_x,domain_y,domain_z] = deal(1:size(field(1).value,1),1:size(field(1).value,2),1:size(field(1).value,3));
        end
end


% Fix axis to domain size and max amplitude (if one isn't already given) for each field
Max=zeros(1,num_fields); Min=Max;
for i=1:num_fields
    
    if field(i).min_amp_given && field(i).max_amp_given
        assert(field(i).max_amp>field(i).min_amp, 'Maximum amplitude must be larger than minimum amplitude')
    end
    
    switch dim
    case 1
        
        % Maximum amplitude for 1D fields
        if field(i).max_amp_given~=0
            Max(i)=field(i).max_amp;
        else
            Max(i)=max(max(field(i).value));
        end
        % minimum amplitude for 1D fields
        if field(i).min_amp_given~=0
            Min(i)=field(i).min_amp;
        else
            Min(i)=min(min(field(i).value));
        end
        assert(~all(all(field(i).value==0)), 'Field "%s" must have values other than 0',field(i).name)
        field(i).axis_size = [domain_x(1) domain_x(end) Min(i) Max(i)];
        
    case 2
        
        % Maximum amplitude for 2D fields
        if field(i).max_amp_given~=0
            Max(i)=field(i).max_amp;
        else
            Max(i)=max(max(max(field(i).value)));
        end
        % Minimum amplitude for 2D fields
        if field(i).min_amp_given~=0
            Min(i)=field(i).min_amp;
        else
            Min(i)=min(min(min(field(i).value)));
        end
        
        % convert domain arrays to matrices if required
        if isvector(domain_x)
            [domain_y,domain_x] = ndgrid(domain_x,domain_y);
        end
        
        if field(i).vector_field
            if strcmp(field(i).plot_type,'quiver3'), domain_z=zeros(size(domain_x)); end
        else
            assert(~all(all(all(field(i).value==0))), 'Field "%s" must have values other than 0',field(i).name)
            field(i).axis_size = [domain_x(1) domain_x(end) domain_y(1) domain_y(end) Min(i) Max(i) Min(i) Max(i)];
        end
        
        
    case 3
        
        % Maximum amplitude for 3D fields
        if field(i).max_amp_given~=0
            Max(i)=field(i).max_amp;
        else
            Max(i)=max(max(max(max(field(i).value))));
        end
        % Minimum amplitude for 2D fields
        if field(i).min_amp_given~=0
            Min(i)=field(i).min_amp;
        else
            Min(i)=min(min(min(min(field(i).value))));
        end
        assert(~all(all(all(all(field(i).value==0)))), 'Field "%s" must have values other than 0',field(i).name)
        field(i).axis_size = [domain_x(1) domain_x(end) domain_y(1) domain_y(end) Min(i) Max(i) Min(i) Max(i)];
        [domain_x,domain_y,domain_z] = ndgrid(domain_x,domain_y,domain_z);
    end
    
    % patch elements have 0 height
    if strcmp(field(i).plot_type, 'patch')
        field(i).patch_colorbar = field(i).axis_size(3:4);
        field(i).axis_size(end) = 0;
    end
end


% Work out suplot parameters
num_vertical_plots = round(sqrt(num_fields));
num_horizontal_plots = round(sqrt(num_fields+num_vertical_plots));


% Set figure up with UI panel
fig = figure;
if fullscreen
    screen_size = get(0,'ScreenSize');
    % fix so that start menu doesn't get in the way
    screen_size(1:2)=screen_size(1:2)+50;  % figure position
    screen_size(3:4)=screen_size(3:4)-120; % figure size
    set(fig,'Position',screen_size);
end
% position = [distance_from_left distance_from_bottom width height]
top_panel = uipanel('Parent', fig,...
    'BackgroundColor','white',...
    'Position',[0 0.1 1 .9]);
bottom_panel = uipanel('Parent', fig,...
    'BackgroundColor','white',...
    'Position',[0 0 1 .1]);


% Create sublot axes
for i=1:num_fields
    ax(i)=subplot(num_vertical_plots,num_horizontal_plots,i,...
        'Parent',top_panel);
end


% Show play putton
handle.play = uicontrol('Parent',bottom_panel,...
    'Style', 'pushbutton',...
    'Position', [20 10 95 20]);


% Show timestep slider
arrow_step = 1/NT;      % steps are percentage of overall range
trough_step = 100/NT;
handle.timeslider = uicontrol('Parent',bottom_panel,...
    'Style', 'slider',...
    'Min',1,'Max',NT,'Value',1,...
    'Position', [120 10 195 20],...
    'Callback', {@select_timestep},...
    'SliderStep',[arrow_step trough_step]);


% Show frame number
k=1;
handle.current_time = uicontrol('Parent',bottom_panel, 'Style','text',...
    'Position',[320 10 45 20],'String',k);


% If save is required
if save_location~=0
    
    % Show 'stop recording' button
    handle.stop_recording = uicontrol('Parent',bottom_panel, 'Style', 'pushbutton',...
        'Position', [370 10 95 20], 'String','Stop Recording');
    
    % Pre-process video file
    writerObj = VideoWriter(save_location);
    
    writerObj.FrameRate = fps;
%     writerObj.FrameRate = 1/frame_delay/2;    % make fps same as previewed in Matlab

    open(writerObj);
    
end


% Play automatically
colormap(customcolor);
play_function



% Play button function
    function play_function(~,~)
        stop_please=0;
        
        % Change play button to stop button
        set(handle.play,'string','Stop',...
            'Callback', @break_function);
        
        % Time loop
        try
            hold on;
            while (k<NT+1 && ~stop_please)
                
                % Subplot all fields given
                for p=1:num_fields
                    plot_field(p)
                end
                
                frame = getframe(fig);
                pause(frame_delay)
                
                % Write a frame to the video file if required
                if save_location~=0
                    set(gca,'nextplot','replacechildren');
                    % set(gcf,'Renderer','zbuffer');     % to work with OpenGL renderer
                    writeVideo(writerObj,frame);
                    
                    % 'stop recording' function
                    set(handle.stop_recording, 'Callback', @stop_recording_fun);
                    
                end
                
                % Show frame number and modify slider position
                set(handle.current_time,'String',k);
                set(handle.timeslider,'Value',k);
                
                % Advance one time step
                k=k+skipTimesteps;
            end
        catch
            warning('Animation window closed before it had stopped')
        end
        
        break_function
    end


% Plot individual field
    function plot_field(i)
        
        % change current axes and clear all previous plots
        axes(ax(i)); cla;
        
        switch field(i).plot_type
            
            case 'plot'
                
                plot(domain_x, field(i).value(:,k));
                axis(field(i).axis_size);
                
            case 'surf'
                
                surf(domain_x,domain_y, field(i).value(:,:,k), 'LineStyle','none');
%                 surf(domain_x,domain_y, field(i).value(:,:,k), 'EdgeColor','none', 'FaceColor','interp');
                view (field(i).dimensions);
                axis(field(i).axis_size);
                if field(i).dimensions==2, colorbar; end
                
            case 'patch'
                
                patches = field(i).patches;
                patches.FaceVertexCData = field(i).value(:,k);
                patch(patches);
                
                view (field(i).dimensions);
                caxis(field(i).patch_colorbar);
                colorbar;
                axis equal;
                
            case 'quiver3'
                
                quiver3(domain_x,domain_y,domain_z, field(i).u(:,:,k),field(i).v(:,:,k),field(i).w(:,:,k));
                view (field(i).dimensions);
                
            case 'quiver'
                
                quiver(domain_x,domain_y, field(i).u(:,:,k),field(i).v(:,:,k));
        end
        
        % plot overlay
        if field(i).plot_overlay
            hold on
            above_animation = field(i).axis_size(end)*ones(length(field(i).overlay),1);
            scatter3(field(i).overlay(:,1),field(i).overlay(:,2),above_animation,10,...
                'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[1 1 1]);
%             plot(field(i).overlay(:,1),field(i).overlay(:,2),...
%                 'MarkerSize',10,'Marker','.','LineStyle','none','Color',[0 0 0])
            hold off
        end
        
        grid(field(i).show_grid)
        title(field(i).name);
    end


% Stop button function
    function break_function(~,~)
        stop_please=1;
        
        % Change stop button to play button
        set(handle.play, 'String', 'Play',...
            'Callback', @play_function);
    end


% Select timestep function
    function select_timestep(hObj,~)
        
        % Stop animation
        break_function
        
        % Get slider value
        k = round(get(hObj,'Value'));
        
        % Subplot all fields given
        for q=1:num_fields
            plot_field(q)
        end
        
        frame = getframe(fig);
        
        % Show frame number
        set(handle.current_time,'String',k);
        set(handle.timeslider,'Value',k);
        
        % Write a frame to the video file if required
        if save_location~=0
            set(gca,'nextplot','replacechildren');
            % set(gcf,'Renderer','zbuffer');     % to work with OpenGL renderer
            writeVideo(writerObj,frame);
            
            % 'stop recording' function
            set(handle.stop_recording, 'Callback', @stop_recording_fun);
            
        end
        
    end


% Stop recording button function
    function stop_recording_fun(~,~)
        stop_please=1;
        
        % Change stop button to play button
        set(handle.play, 'String', 'Play',...
            'Callback', @play_function);
        
        % Fade stop recroding button
        set(handle.stop_recording, 'Enable','off');
        
        % Close the video file once video has finished
        close(writerObj);
        
        % alternate method to create a movie:
        % movie2avi(M,save_location, 'compression', 'None');
        
        save_location=0;
    end


end
