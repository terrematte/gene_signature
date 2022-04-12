library("tidyverse")

pubmed <- read_csv("data/pubMed_timeline_results_by_year.csv", skip = 1)
pubmed <- as.data.frame(pubmed)
pubmed <- pubmed[pubmed$Year < 2021 ,]
sum(pubmed[pubmed$Year < 2021 ,"Count"])

sum(pubmed[pubmed$Year < 2021 & pubmed$Year >= 2015,"Count"])

sum(pubmed[pubmed$Year < 2021 & pubmed$Year >= 2015,"Count"])


pubmed <- rbind(pubmed, data.frame(Year=setdiff(c(2008:2020), pubmed$Year), Count = 0))

pubmed$inclusion <- "Excluded"
pubmed$inclusion[pubmed$Year < 2021 & pubmed$Year >= 2015] <- "Included"

pdf("figs/figA01_pubmed_search.pdf", width = 10, height = 5)
ggplot(data=pubmed, aes(x=Year, y=Count))  +
  #geom_line(color="#69b3a2", size=2) +
  #geom_line(color="#69b3a2", size=2) +
  geom_line(mapping = aes(x = ifelse(Year<=2015 , Year, 0)), color = "pink",  size=2) +
  geom_line(mapping = aes(x = ifelse(Year>=2015 & Year < 2021 , Year, 0)), color = "#69b3a2", size=2) +
  geom_label(label=pubmed$Count, 
             nudge_x =0, nudge_y = 2) + 
  geom_area(mapping = aes(x = ifelse(Year>=2015 & Year < 2021 , Year, 0)), fill = "#69b3a2", alpha=0.4) +
  geom_area(mapping = aes(x = ifelse(Year<=2015 , Year, 0)), fill = "pink",  alpha=0.4) +
  geom_point() +
  #scale_x_continuous(breaks=custom_breaks) + +
  #scale_x_continuous(breaks = c(2008:2020), labels = c(2008:2020)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45)) +
  #scale_x_continuous(breaks = integer_breaks()) 
  #scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) +
  #scale_y_continuous(breaks=c(0,2,4,6,8,10,12)) +
  scale_x_continuous(breaks=rev(pubmed$Year), limits = c(2008,2020)) +
  ylab("Results") 
dev.off()
