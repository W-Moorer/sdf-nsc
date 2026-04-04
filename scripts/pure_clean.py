import sys
import codecs
with open('papers/sections/03_methodology.tex', 'rb') as f:
    b = f.read()

# decode safely
text = b.decode('utf-8', errors='ignore')

# remove null bytes and garbled sections
text = text.replace('\0', '')

# identify the garbled text and remove it
text_to_replace = 'subsectionWNbQb gR:SNAmbM v(uirteHh(T opology −AwareContactM a）'
text_to_replace2 = 'subsectionWNbQb gR:SNAmbM v(uirteHh(T opology -AwareContactM a）'
text_to_replace3 = 'subsectionWNbQb gR:SNAmbM v(uirteHh(T opology -AwareContactM a'

if text_to_replace in text: print("Found garbled 1")
elif text_to_replace2 in text: print("Found garbled 2")
elif text_to_replace3 in text: print("Found garbled 3")

import re
text = re.sub(r'\\subsectionWNbQb\s*.*?a）', '', text)
text = re.sub(r'subsectionWNbQb\s*.*?a）', '', text)
text = re.sub(r'subsectionWNbQb\s*.*?a', '', text)
text = re.sub(r'WNbQb gR:SNAmbM v\(uirteHh\(T opology -AwareContactM a', '', text)
text = re.sub(r'WNbQb gR:SNAmbM v\(uirteHh\(T opology −AwareContactM a', '', text)
text = re.sub(r'WNbQb.*AwareContactM a[）\)]*', '', text)


# Now, search for the end of 3.4
# We know the garbled text is at the end of 3.4.
# We also know "Scheme A" "方案A" needs to be removed.
# 3.4 is \subsection{基于拓扑结构分区与流形降维的通用物理计算方案 (Topology-Aware Contact Manifold Subsumption)}
# Let's remove the whole section or modify it 

# write it back cleanly
with open('papers/sections/03_methodology.tex', 'w', encoding='utf-8') as f:
    f.write(text)

with open('clean_test.txt', 'w', encoding='utf-8') as f:
    f.write(text[-1500:])

