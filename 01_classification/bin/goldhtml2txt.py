#!/usr/bin/python3

import sys

tablefound=False
cellfound=False
values=[]
for line in sys.stdin:

  if not tablefound:
    if line.strip() == '<tr class="odd">':
      tablefound=True

  else:
    if line.strip() == '</tbody>':
      tablefound=False
    else:
      if not cellfound:
        if line.strip() == '<td>':
          cellfound=True
          values.append('')
        elif line.strip() == '</tr>':
          sys.stdout.write("%s\n" % "\t".join(values))
          values=[]
      else:
        if line.strip() == '</td>':
          cellfound=False
        else:
          values[-1] += line.strip()
