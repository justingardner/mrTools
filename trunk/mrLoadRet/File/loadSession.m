function [session, groups, version] = loadSession(dirPathStr)
% function [session, groups, mrLoadRetVersion] = loadSession([dirPathStr])
%
% Loads the mrSESSION and groups structures from mrSESSION.mat.
%
% dirPathStr: directory that contains the mrSESSION.mat file
% defaults to current directory (pwd)
%
% djh, 2/17/2001
% 7/16/02 djh, update mrSESSION to version 3.01
% 11/2004 djh, update to version 4.0

% current version of mrLoadRet
version = mrLoadRetVersion;

% Default dirPathStr = pwd
curDir = pwd;
if ieNotDefined('dirPathStr')
    dirPathStr = curDir;
end
pathStr = fullfile(dirPathStr,'mrSession.mat');

% Load mrSESSION.mat to get session and groups
if exist(pathStr,'file')
    load(pathStr)
    % Check version numbers for consistency
    if (session.mrLoadRetVersion < version)
        mrErrorDlg(['This session was created with mrLoadRet version ',...
            num2str(session.mrLoadRetVersion),...
            '. Current mrLoadRet version is ',...
            num2str(version),...
            '. Use mrInitRet and start from scratch.']);
    end
else
%    mrErrorDlg(['No mrSESSION.mat file in directory: ',dirPathStr]);
    session = [];
    groups = [];
end

