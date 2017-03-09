function [nx,ny] = boxshape(x,y,xmin,xmax,ymin,ymax)
    if nargin < 6
        ymax = 1;
    end
    if nargin < 5
        ymin = 0;
    end
    if nargin < 4
        xmax = 1;
    end
    if nargin < 3
        xmin = 0;
    end

    if size(x,2)==1
        % x = [x;x(1)];
        % y = [y;y(1)];
    %     nxy = [max(xmin,min(xmax,x)),max(ymin,min(ymax,y))];
    %     nxy = [x,y];
    %     nxy = NaN(numel(x)*5,2);
        cdx = circdiff(x);
        cdy = circdiff(y);
        nmatch = abs(cdx)>0;

        rx = repmat(x',5,numel(x));
        ry = repmat(y',5,numel(y));
        nxy = [rx(:),ry(:)];
    %     nxy(1:5:end,1) = x;
    %     nxy(1:5:end,2) = y;
    %     nxy(2:2:end-1,1) = x(2:end);
    %     nxy(2:2:end-1,2) = y(2:end);
    %     nxy(end,:) = [x(1),y(1)];
        xpost = [x(2:end);x(1)];
        ypost = [y(2:end);y(1)];
    %     pre = [x(end),y(end);x(1:end-1),y(1:end-1)];
    %     xrev = xcur < xpost;
    %     yrev = ycur < ypost;
    %     xcur(xrev) = xpost(xrev);
    %     ycur(yrev) = ypost(yrev);
    %     xpost(xrev) = x(xrev);
    %     ypost(yrev) = y(yrev);
    %     m = (ypost-y)./(xpost-x);
    %     c = (y-ypost)./(m.*x);
        cur = [x,y];
        post = [xpost,ypost];

        for i = 1:numel(x)
            c = 1+(i-1)*5;
            crossval(xmin,false,1);
            crossval(xmax,false,2);
            crossval(ymin,true,3);
            crossval(ymax,true,4);
        end
    %     c = 1;
    %     i = 1;
    %     while c <= size(nxy,2)
        %     if x(i) < xmin || x(i) > xmax
        %         if x(i) < xmin
        %             nx(c) = xmin;
        %         else
        %             nx(c) = xmax;
        %         end
        %         
        %         ny(c) = ypre(c)+g(i)*(nx(c)-xpre(i));
        %     end

    %     crossval(xmin,false,false);
    %     crossval(xmax,false,false);
    %     crossval(ymin,true,false);
    %     crossval(ymax,true,true);

        %     if ny(i) < ymin || ny(i) > ymax
        %         if y(i) < ymin
        %             ny(i) = ymin;
        %         else
        %             ny(i) = ymax;
        %         end
        %         
        %         
        %         nx(c) = xpre(i)+(ny(c)-ypre(i))/g(i);
        %     end

    %         i = i+1;
    %         c = c+1;
    %     end

    %     nx = min(xmax,max(xmin,nxy(:,1)));
    %     ny = min(ymax,max(ymin,nxy(:,2)));
        nx = nxy(:,1);
        ny = nxy(:,2);
    %     nx = nx(~isnan(nx));
    %     ny = ny(~isnan(ny));
        torem = isnan(nx) | isnan(ny) | nx < xmin | nx > xmax | ny < ymin | ny > ymax;
    %     nx = min(xmax,max(xmin,nx(~isnan(nx))));
    %     ny = min(ymax,max(ymin,ny(~isnan(ny))));
        [nx,ny] = shapedeldup(nx(~torem),ny(~torem));
%         norep = (circdiff(nx)~=0 | circdiff(ny)~=0);
%         nx = nx(norep);
%         ny = ny(norep);
    else
        cx = cell(1,size(x,2));
        cy = cell(1,size(x,2));
        len = -inf;
        for i = 1:size(x,2)
            [cx{i},cy{i}] = boxshape(x(:,i),y(:,i),xmin,xmax,ymin,ymax);
            if length(cx{i}) > len
                len = length(cx{i});
            end
        end
%         [nx,ny] = deal(NaN(len,size(x,2)));
%         for i = 1:size(x,2)
%             ind = 1:length(cx{i});
%             nx(ind,i) = cx{i};
%             ny(ind,i) = cy{i};
%             ind2 = length(cx{i})+1:len;
%             nx(ind2,i) = nx(length(cx{i}),i);
%             ny(ind2,i) = ny(length(cy{i}),i);
%         end
        [nx,ny] = appendshapes(cx{:},cy{:});
%         keyboard
    end

    function crossval(cval,isy,ioff)
        xy = isy+1;
        xy2 = 2-isy;
%         if dobreak
%             fprintf('lastx: %f\nlasty: %f\n',pre(i,:));
%             fprintf('x: %f\ny: %f\n\n',cur(i,:));
%         end
        pos = (post(i,xy) >= cval && cur(i,xy) < cval);
        neg = (post(i,xy) < cval && cur(i,xy) >= cval);
        if pos || neg
%             cxy = [NaN NaN];
%             cxy(xy) = cval;
%             cxy(xy2) = 
            cxy = cval;
            chg = (cval-cur(i,xy))./(post(i,xy)-cur(i,xy));
            cxy2 = cur(i,xy2)+chg*(post(i,xy2)-cur(i,xy2));
            nxy(c+(ioff:4),xy) = cval;
%             if isy
%                 grad = m(i);
%             else
%                 grad = 1/m(i);
%             end
            nxy(c+(ioff:4),xy2) = cxy2;
%             nxy = [nxy(1:c+ioff,:);cxy;nxy(c+1+ioff:end,:)];
%             c = c+1;
%             disp(1)
        end
    end
end