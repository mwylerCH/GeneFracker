pdf(NULL)
args <- commandArgs(TRUE)
TABLE <- args[1]
OUTPLOT <- args[2]
suppressMessages(suppressWarnings(require(tidyverse)))

tabella <- read.delim(TABLE, sep = '\t', header = T)

tabella %>%
  mutate(MidPoint = (Start+End)/2,
        # geneName = factor(geneName, levels = unique(geneName)),
         geneLength = geneLength) %>% 
  ggplot()+
  geom_segment(aes(y=geneName, yend = geneName,
               x=1 , xend = geneLength))+
  geom_rect(aes(xmin = Start, xmax = End,
                ymin = as.numeric(geneName)-0.2, ymax = as.numeric(geneName)+0.2,
                fill = Domain))+
  geom_text(aes(y=geneName, x = MidPoint, label = Domain), size = 2)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = 'black'))+
  labs(y = '',
       x = '')
suppressMessages(suppressWarnings(ggsave(OUTPLOT, width = 10)))
