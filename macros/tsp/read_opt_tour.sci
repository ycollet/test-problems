function OptTour = read_opt_tour(Filename)
[fd, err] = mopen(Filename,'r');
if (err) then
  printf('error while reading file %s\n',Filename);
  abort;
end

Line = mgetl(fd,1);
Line = stripblanks(tokens(Line,':'));

while(Line(1)~='TOUR_SECTION')
  Line = mgetl(fd,1);
  Line = stripblanks(tokens(Line,':'));
end

Line = mgetl(fd,1);
Line = stripblanks(Line);
Index = 0;
while(Line(1)~='EOF')
  Index = Index + 1;
  OptTour(Index) = eval(Line(1));
  Line = mgetl(fd,1);
  Line = stripblanks(Line);
end // While
OptTour = OptTour(1:$-1); // Because the last index is -1 to tell to close the loop
mclose(fd);
endfunction
