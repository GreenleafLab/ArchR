---
name: Bug/Error Report - NO USAGE QUESTIONS / FEATURE REQUESTS (Use Discussions!)
about: Create a bug/error report to help us improve ArchR. NOT to be used for usage
  questions or feature requests!
title: ''
labels: bug
assignees: ''

---

This is an issue template made by the developers of ArchR. You MUST follow these instructions.

Questions related to how to use ArchR or requests for new features should be posted in the Discussions forum (https://github.com/GreenleafLab/ArchR/discussions).

Before you submit this Bug Report please update ArchR to the latest stable version and make sure that this issue has not already been fixed in the latest release. ArchR is still in active development and we will fix problems as they arise. To update ArchR:

devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

If your issue persists, then please submit this bug report.

PLEASE FILL OUT THE RELEVANT INFORMATION AND DELETE THE UNUSED PORTIONS OF THIS ISSUE TEMPLATE.

**Attach your log file**
ArchR has a built-in logging functionality for all complex functions. You MUST attach your log file (indicated in the console output) to this issue. Just drag and drop it here.

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
To help us optimally address your issue, please try to reproduce this issue using the tutorial hematopoiesis dataset and provide us the command(s) to reproduce your bug. Our first question to you will be "can you reproduce this with the tutorial dataset" so please do this.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem. Do not screenshot code or text but embed this in markdown using triple-backticks.

**Session Info**
If you do not have a log file because the function that caused the error does not produce one, please paste the output of "sessionInfo()" here.

**Additional context**
Add any other context about the problem here.
