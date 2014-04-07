function [insideNodes,insideNodeStruct] = floodFillFindPerim(mesh,perimDist,startNode,statusHandle)
% A flood-fill method to find the mesh perimeter
%
% [insideNodes,insideNodeStruct]=floodFillFindPerim(mesh,perimDist,startNode,busyHandle)
%
% Find expanding rings connected to the start node. 
% Stop when all the members of the ring exceed perimDist from the startNode
% But cleverer than this :)
%
% Also stores the average distance of each new set and its offset - May be
% useful for tacking down points later This can break quite easily - for
% example when the floodfill runs around all sides of a bump. But we can
% fix this later  ...
%
% AUTHOR:  Wade
% DATE : 020701 last modified
 
compute=1;
recompute=0;

while compute
  withinDist = 1;
  insideNodes   = startNode;
  currentNodes  = startNode;
  counter       = 0;
  insideNodeStruct.offset = 0;
  insideNodeStruct.avDist = 0;
  nVerts = length(mesh.connectionMatrix);
  cmBackup = mesh.connectionMatrix;
  % Limit to generate no more that 10000 rings....
  while ((withinDist) && (counter<10000)) 

     % What nodes are connected to the current ones?
     [newRows connected]=find(mesh.connectionMatrix(currentNodes,:));

     if (~isempty(connected))

        connected=unique(connected(:));   
        insideNodes=[insideNodes;connected];

        currentNodes=connected;
        notCnodes=setdiff(1:nVerts,currentNodes);

        % This is a cute way of zeroing columns in a sparse matrix. Much
        % faster than foo(:,currentNodes)=0;  
        diagMat=sparse(notCnodes,notCnodes,ones(length(notCnodes),1),nVerts,nVerts);
        mesh.connectionMatrix=(mesh.connectionMatrix)*(diagMat);

        nodeDists=mesh.dist(currentNodes);
        if recompute  % if the flood fill has gone too far the first time, we change the criterion 
                      % to stop as soon as the average perimeter is reached  
          withinDist=mean(nodeDists(:))<perimDist;
        else
          withinDist=any(nodeDists<perimDist);
        end

        % Average distance
        insideNodeStruct.avDist=[insideNodeStruct.avDist;mean(nodeDists(:))];

        % Offset to the list of current rings
        insideNodeStruct.offset=[insideNodeStruct.offset;length(connected)];
  %       round([nnz(mesh.connectionMatrix) insideNodeStruct.avDist(counter+2) insideNodeStruct.offset(counter+2) nnz(nodeDists<perimDist)])

     else
           withinDist=0;
     end   

     counter=counter+1; 
  end
  if recompute
    str=sprintf('** Warning: new max average perimeter:%.2f mm\n',insideNodeStruct.avDist(end));
    statusStringAdd(statusHandle,str);
    compute=0;
  else
    if insideNodeStruct.avDist(end)>perimDist %if the last ring found is more than the perimeter asked
      beep;
      str=sprintf('** Warning: flood fill went beyond perimeter (max average: %.2f mm). Have you checked/corrected for topological defects?',insideNodeStruct.avDist(end));
      statusStringAdd(statusHandle,str);
      if insideNodeStruct.avDist(end)>perimDist*1.05 %if the last ring found is more than 105% the perimeter asked
        recompute=1;
        mesh.connectionMatrix = cmBackup;
        statusStringAdd(statusHandle,sprintf('** Warning: Recomputing with different stopping criterion...'));
      else
        compute=0;
      end
    else
      compute=0;
    end
  end
end

insideNodes=unique(insideNodes);

return;
