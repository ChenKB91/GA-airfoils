function varargout = parsec_gui(varargin)
% PARSEC_GUI MATLAB code for parsec_gui.fig
%      PARSEC_GUI, by itself, creates a new PARSEC_GUI or raises the existing
%      singleton*.
%
%      H = PARSEC_GUI returns the handle to a new PARSEC_GUI or the handle to
%      the existing singleton*.
%
%      PARSEC_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARSEC_GUI.M with the given input arguments.
%
%      PARSEC_GUI('Property','Value',...) creates a new PARSEC_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before parsec_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to parsec_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help parsec_gui

% Last Modified by GUIDE v2.5 27-Feb-2019 11:17:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @parsec_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @parsec_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before parsec_gui is made visible.
function parsec_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to parsec_gui (see VARARGIN)

% Choose default command line output for parsec_gui
% handles.preset=[0.015 0.45 0.07 -1.75 0.3 -0.06 0.45 0 0 13.63 32.16];

handles.preset=zeros(6,11);
handles.preset(1,:) = [0.039134 0.333752 0.150554 -0.888058 0.545012 -0.065095 2.399356 -0.071929 0.000000 25.001556 -4.880719];
% handles.preset(1,:) = [0.037402 0.333752 0.009094 -0.187713 0.572379 -0.072605 1.017089 -0.065221 0.000000 20.325262 -1.184745];
handles.preset(2,:) = [0.032224 0.465039 0.059102 -2.303902 0.321705 -0.058355 1.126030 -0.009451 0.000000 24.601266 -5.791860];
handles.preset(3,:) = [0.032224 0.465039 0.059102 -0.963002 0.559825 -0.058355 1.700247 -0.009451 0.000000 24.601266 3.522452];
handles.preset(4,:) = [0.104000 0.439291 0.043900 -2.518931 0.367363 -0.065095 0.424392 -0.071929 0.000000 26.801602 -3.566474];
handles.preset(5,:) = [0.037402 0.333752 0.009094 -0.888058 0.420847 -0.072605 0.435196 -0.071929 0.000000 26.801602 -5.891533];
handles.preset(6,:) = [0.037402 0.415993 0.009094 -0.888058 0.559825 -0.072605 0.424392 -0.071929 0.000000 22.868813 4.177060];
handles.preset(7,:) = [0.037402 0.333752 0.071370 -1.074357 0.572379 -0.065095 1.704784 -0.071929 0.000000 22.150314 3.522452];
handles.preset(8,:) = [0.037402 0.333752 0.009094 -0.547398 0.491267 -0.072605 1.416645 -0.071929 0.000000 22.868813 -5.891533];
% handles.preset=[0.01 0.3 0.06 -0.45 0.3 -0.06 0.45 0 0 0 17];
% handles.preset=[0.0146,0.3025,0.06,-0.4928,0.3016,0.06,-0.4848,-0.1039,0,-2.7791,9.2496];
% handles.variance=[
handles.output = hObject;

handles.min = [0.024793 0.25 0.07 -0.8244 0.26 -0.08 0.4244 -0.05 0.015 2 -60];
handles.max = [0.024793 0.35 0.08 -0.4244 0.34 -0.07 0.8244 0.05 0.015 18 60];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes parsec_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = parsec_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in plotbutton.
function plotbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plotbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    rle = str2double(get(handles.rle,'string'));
    xup = str2double(get(handles.xup,'string'));
    yup = str2double(get(handles.yup,'string'));
    yxxup = str2double(get(handles.yxxup,'string'));
    xlo = str2double(get(handles.xlow,'string'));
    ylo = str2double(get(handles.ylow,'string'));
    yxxlow = str2double(get(handles.yxxlow,'string'));
    yte = str2double(get(handles.yte,'string'));
    deltayte = str2double(get(handles.deltayte,'string'));
    alpha = str2double(get(handles.alpha,'string'));
    beta = str2double(get(handles.beta,'string'));
    
    pltoptions = ["-","o","o-","."];
    chosen_option = pltoptions(get(handles.popupmenu1,'value'));
% cla(handles.axes1)
%     disp(pltoptions(get(handles.popupmenu1,'value')))

    par = [rle,xup,yup,yxxup,xlo,ylo,yxxlow,yte,deltayte,alpha,beta];
    pt_distribute = get(handles.popup_distribution,'value');
    switch pt_distribute
        case 1
            [pts,self_cross] = evenpar(par);
%             disp(pts)
            plot(handles.axes1,pts(:,1),pts(:,2),chosen_option)
            hold on
            polyin = polyshape({pts(:,1)},{pts(:,2)});
            [comx,comy] = centroid(polyin);
        case 2
            [pts,self_cross] = parsecpoints(par);
%             disp(pts)
            plot(handles.axes1,pts.x,pts.y,chosen_option)
            hold on
            polyin = polyshape({pts.x},{pts.y});
            [comx,comy] = centroid(polyin);
    end
    if self_cross
        set( handles.selfcrosstext,'String','Warning: Self-Crossing!');
    else
        set( handles.selfcrosstext,'String','');
    end
%      disp([comx comy])
     
    if get(handles.buttonPlotCOM,'value')
        plot(handles.axes1,comx,comy,'r*')
        hold on
    end
    
    axis equal

function xup_Callback(hObject, eventdata, handles)
% hObject    handle to xup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xup as text
%        str2double(get(hObject,'String')) returns contents of xup as a double


% --- Executes during object creation, after setting all properties.
function xup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yup_Callback(hObject, eventdata, handles)
% hObject    handle to yup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yup as text
%        str2double(get(hObject,'String')) returns contents of yup as a double


% --- Executes during object creation, after setting all properties.
function yup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yxxup_Callback(hObject, eventdata, handles)
% hObject    handle to yxxup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yxxup as text
%        str2double(get(hObject,'String')) returns contents of yxxup as a double


% --- Executes during object creation, after setting all properties.
function yxxup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yxxup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xlow_Callback(hObject, eventdata, handles)
% hObject    handle to xlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlow as text
%        str2double(get(hObject,'String')) returns contents of xlow as a double


% --- Executes during object creation, after setting all properties.
function xlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rle_Callback(hObject, eventdata, handles)
% hObject    handle to rle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rle as text
%        str2double(get(hObject,'String')) returns contents of rle as a double


% --- Executes during object creation, after setting all properties.
function rle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ylow_Callback(hObject, eventdata, handles)
% hObject    handle to ylow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ylow as text
%        str2double(get(hObject,'String')) returns contents of ylow as a double


% --- Executes during object creation, after setting all properties.
function ylow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ylow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yxxlow_Callback(hObject, eventdata, handles)
% hObject    handle to yxxlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yxxlow as text
%        str2double(get(hObject,'String')) returns contents of yxxlow as a double


% --- Executes during object creation, after setting all properties.
function yxxlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yxxlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yte_Callback(hObject, eventdata, handles)
% hObject    handle to yte (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yte as text
%        str2double(get(hObject,'String')) returns contents of yte as a double


% --- Executes during object creation, after setting all properties.
function yte_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yte (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function deltayte_Callback(hObject, eventdata, handles)
% hObject    handle to deltayte (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deltayte as text
%        str2double(get(hObject,'String')) returns contents of deltayte as a double


% --- Executes during object creation, after setting all properties.
function deltayte_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deltayte (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_Callback(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta as text
%        str2double(get(hObject,'String')) returns contents of beta as a double


% --- Executes during object creation, after setting all properties.
function beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FillPreset.
function FillPreset_Callback(hObject, eventdata, handles)
% hObject    handle to FillPreset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

textareas = [...
handles.rle
handles.xup
handles.yup
handles.yxxup
handles.xlow
handles.ylow
handles.yxxlow
handles.yte
handles.deltayte
handles.alpha
handles.beta];
dat = str2double(split(get(handles.edit13,'string')));
for i = 1:11
    set(textareas(i),'string',dat(i))
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_distribution.
function popup_distribution_Callback(hObject, eventdata, handles)
% hObject    handle to popup_distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_distribution contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_distribution


% --- Executes during object creation, after setting all properties.
function popup_distribution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exitbutton.
function exitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exit= questdlg('Confirm Exit?', 'Exit', 'Yes','No','No'); 
if exit =='Yes'
    close(gcf ) 
end


% --- Executes on button press in FillRandom.
function FillRandom_Callback(hObject, eventdata, handles)
% hObject    handle to FillRandom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rng('shuffle')
min = handles.min; % [0   0.3 0   -3 0.2 -0.3 0  -0.1 0 -10 20];
max = handles.max;% [0.1 0.6 0.3 0   0.6  0  3  0.1 0  10 30];
% min = [0.02  0.32 0.077 -0.65 0.15 -0.05 -0.60 0 0.01 -4.55 15];
% max = [0.023 0.37 0.080 -0.63 0.19 -0.02 -0.75 0 0.01 -4.90 15.1];
r = min + (max-min).*rand(1,11);
textareas = [...
handles.rle
handles.xup
handles.yup
handles.yxxup
handles.xlow
handles.ylow
handles.yxxlow
handles.yte
handles.deltayte
handles.alpha
handles.beta];
for i = 1:11
    set(textareas(i),'string',r(i))
end
plotbutton_Callback(handles.plotbutton, eventdata, handles)


% --- Executes on button press in buttonPlotCOM.
function buttonPlotCOM_Callback(hObject, eventdata, handles)
% hObject    handle to buttonPlotCOM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttonPlotCOM



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hold on

% --- Executes on button press in offbutton.
function offbutton_Callback(hObject, eventdata, handles)
% hObject    handle to offbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1)