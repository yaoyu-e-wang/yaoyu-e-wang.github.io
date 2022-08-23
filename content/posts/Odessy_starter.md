---
title: "Odessy_starter"
date: 2020-06-24T15:08:42-04:00
draft: false
---
This note provides introduction guidance to the HUIT Odessy system with FAQ 


## How do I login to Odyssey

ssh into yewang@login.rc.fas.harvard.edu with password and with OpenAuth DUO set up.  The information is at https://www.rc.fas.harvard.edu/resources/documentation/openauth/


## Center Shared Directory

The center shared directory is located at:

/n/qbrc_center 

and can be mounted (https://www.rc.fas.harvard.edu/resources/documentation/mounting-storage/) as smb://qbrc.rc.fas.harvard.edu/qbrc_center

## Update members within qbrc_center resources

If users with preexisting RC accounts
want to join your group, they can do that with Portal:

https://portal.rc.fas.harvard.edu/request/grants/add

Otherwise, new users can list you in the PI section when signing up for an
account. Both will send you an email notification to approve/reject their
request through Portal.

de


## Transfer Data from external source to fas

Here is the link to fas FAQ page:

https://www.rc.fas.harvard.edu/resources/documentation/linux/rsync/

Use example:

rsync -avz  foo/ MYUSERNAME@odyssey.fas.harvard.edu:~/foo/







Since your storage space and Fairshare (https://www.rc.fas.harvard.edu/fairshare/)
is part of the qbrc_center lab/unix group, anyone who wants to use those
resources will need to be in that group. If users with preexisting RC accounts
want to join your group, they can do that with Portal (https://portal.rc.fas.harvard.edu/request/grants/add).
Otherwise, new users can list you in the PI section when signing up for an
account. Both will send you an email notification to approve/reject their
request through Portal.

If you or your lab/staff have any needs or issues, please have a look at our
documentation (https://rc.fas.harvard.edu). We can also be contacted via email (http://rchelp@rc.fas.harvard.edu)
or tickets (http://portal.rc.fas.harvard.edu/).


