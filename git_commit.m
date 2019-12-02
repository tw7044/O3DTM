% GIT_COMMIT creates Git commit and pushes changes to github

commit_msg = input('Commit message: ', 's');

disp('> cd ..');
cd ..

disp('> !git add .')
!git add .

disp(strcat('> !git commit -m "', commit_msg, '"'));
system(strcat('git commit -m "', commit_msg, '"'));

disp('> !git push github --all');
!git push github --all

disp('> !git push github --tags');
!git push github --tags

disp('> cd code');
cd code