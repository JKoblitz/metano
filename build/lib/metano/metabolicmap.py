from __future__ import print_function
from builtins import str
from builtins import map
from builtins import object
import sys
from itertools import chain
from pprint import pprint

import MySQLdb
from _mysql_exceptions import OperationalError
import networkx as nx

def flatten(l):
    return list(chain.from_iterable(l))

def listToString(l):
    return "', ".join("'" + str(i) for i in l) + "'"

class MetabolicMap(object):

    reactionShape = "ellipse"
    reactionColor = "#ffcd00"
    
    substanceShape = "roundrectangle"
    substanceColor = "#89a400"
    
    edgeWidth = "3"
    edgeColor =	"#000000"
    edgeArrow = "standard"
    
    fontSize = "14"

    # TODO faded colors    
#    reactionColor = "#ffcd0030"
#    substanceColor = "#89a40030"

    def __init__(self, server, user, password, database, verbosityLevel=0):
        
        self.server = server
        self.user = user
        self.__password = password
        self.database = database
        self.verbosityLevel = verbosityLevel
        self.__connection = None
        self.clearGraph()
    
    def __del__(self):
        self.closeConnection()
    
    def __getConnection(self):
        
        if self.__connection is None: # TODO check for timeout
            # open new db connection
            self.__connection = MySQLdb.connect(self.server, self.user,
                                                self.__password, self.database)
        return self.__connection

    def closeConnection(self):
        
        if self.__connection is not None:
            try:
                # close db connection
                self.__connection.commit()
                self.__connection.close()
                del self.__connection
                self.__connection = None
            except OperationalError:
                # connection is already closed
                del self.__connection
                self.__connection = None

    def clearGraph(self):
        self.graph=nx.DiGraph()

    def printStats(self):
            print("graph has",self.graph.number_of_nodes(),"nodes and", \
                  self.graph.number_of_edges(),"edges")    

    # TODO try other export formats
    def writeToGML(self, filename):
        if self.verbosityLevel > 0:
            self.printStats()
        nx.write_gml(self.graph, filename)
        
#    def writeToDot(self, filename):
#        if self.verbosityLevel > 0:
#            self.printStats()
#        nx.write_dot(self.graph, filename)

#    def writeToGraphML(self, filename):
#        # TODO fails if there are no edges
#        if self.verbosityLevel > 0:
#            self.printStats()
#        nx.write_graphml(self.graph, filename)

    def drawReactionNodes(self, reactionNodeList):

        connection = self.__getConnection()
        cursor = connection.cursor()

        reactionNodeString = listToString(reactionNodeList)
        
        # add all edges and nodes to the graph
        
        cursor.execute("""
        SELECT DISTINCT edge.source_node_id, edge.target_node_id
        FROM edge
        WHERE (edge.source_node_id IN (%s)
        OR edge.target_node_id IN (%s))
        """ % (reactionNodeString, reactionNodeString))
        
        result = cursor.fetchall()
        
        if len(result)<=0:
            print("No edges to draw.")
            return

        self.graph.add_edges_from(result)

        # assign edge attributes
        for start, end  in self.graph.edges_iter():

            # skip already existing edges
            if "graphics" in self.graph[start][end]:
                continue
            
            self.graph[start][end]["graphics"]= \
                { "width":float(self.edgeWidth),
                 "fill":self.edgeColor,
                 "targetArrow":self.edgeArrow
                 }

        # list of all node ids as string
        nodeString = listToString(flatten(result))
        
        # get position, width and height of each node
        cursor.execute("""
        SELECT node_id, x_pos, y_pos
        FROM node
        WHERE node_id IN (%s)
        """ % nodeString)

        result = cursor.fetchall()
        positionDict = dict( (nodeID, (x, y)) for nodeID, x, y in result)

        # get EC numbers as labels
        cursor.execute("""
        SELECT enzyme2node.node_id, enzyme.ec_number
        FROM enzyme2node INNER JOIN enzyme ON (enzyme.enzyme_id=enzyme2node.enzyme_id)
        WHERE enzyme2node.node_id IN (%s)
        """ % nodeString)
        result = cursor.fetchall()
        
        # assign reaction node attributes
        for nodeID, label in result:
            
            # skip already existing nodes
            if "graphics" in self.graph.node[nodeID]:
                continue
            
            self.graph.node[nodeID]["label"]=label
            self.graph.node[nodeID]["graphics"]= \
                {"x":positionDict[nodeID][0],
                 "y":positionDict[nodeID][1],
                 "w":max(30.0, len(label)*9.), # TODO magic numbers
                 "h":30.0,
                 "type":self.reactionShape,
                 "fill":self.reactionColor
                 }
            self.graph.node[nodeID]["LabelGraphics"] = \
                {"fontSize":int(self.fontSize)}

        # get substance names as labels
        cursor.execute("""
        SELECT substance2node.node_id, substance.name
        FROM substance2node INNER JOIN substance ON (substance.substance_id=substance2node.substance_id)
        WHERE substance2node.node_id IN (%s)
        """ % nodeString)
        result = cursor.fetchall()

        # assign substance nodes attributes        
        for nodeID, label in result:
            
            # skip already existing nodes
            if "graphics" in self.graph.node[nodeID]:
                continue
            
            self.graph.node[nodeID]["label"]=label
            self.graph.node[nodeID]["graphics"]= \
                {"x":positionDict[nodeID][0],
                 "y":positionDict[nodeID][1],
                 "w":max(30.0, len(label)*9.), # TODO magic numbers
                 "h":30.0,
                 "type":self.substanceShape,
                 "fill":self.substanceColor
                 }
            self.graph.node[nodeID]["LabelGraphics"] = \
                {"fontSize":int(self.fontSize)}
        

    # TODO MetaCyc only
    def drawReactions(self, reactionIDlist):
        
        connection = self.__getConnection()
        cursor = connection.cursor()
        
        # translate reaction IDs to enzyme nodes
        
        reactionString = listToString(reactionIDlist)
        
        cursor.execute("""
        SELECT r2n.node_id
        FROM reaction2kegg_metacyc_brenda AS r2db
        INNER JOIN reaction2enzymenode AS r2n
        ON (r2db.reaction_id=r2n.reaction_id)
        WHERE Reaction_ID_MetaCyc IN (%s)
        """ % reactionString)
        
        result = cursor.fetchall()

        self.drawReactionNodes(list(map(str,flatten(result))))
    
    def drawEnzymes(self, ecList):
        
        connection = self.__getConnection()
        cursor = connection.cursor()
        
        ecString = listToString(ecList)
        
        # translate EC numbers to enzyme nodes
        
        cursor.execute("""
        SELECT enzyme2node.node_id
        FROM enzyme INNER JOIN enzyme2node ON (enzyme.enzyme_id=enzyme2node.enzyme_id)
        WHERE enzyme.ec_number IN (%s)
        """ % ecString)
        
        result = cursor.fetchall()
        
        self.drawReactionNodes(list(map(str,flatten(result))))
        
        # TODO return list of drawn enzymes
    
    def drawReactionEnzyme(self, reactionEClist):

        connection = self.__getConnection()
        cursor = connection.cursor()
            
        reactionString = listToString(rid for rid, _ in reactionEClist)
        ecString = listToString(ec for _, ec in reactionEClist if ec is not None)
        
        cursor.execute("""
        -- translate reaction id to node id
        SELECT r2n.node_id
        FROM reaction2kegg_metacyc_brenda AS r2db
        INNER JOIN reaction2enzymenode AS r2n
        ON (r2db.reaction_id=r2n.reaction_id)

        -- translate node id to ec numbers
        INNER JOIN enzyme2node AS e2n
        ON (r2n.node_id=e2n.node_id)
        INNER JOIN enzyme AS e
        ON (e2n.enzyme_id=e.enzyme_id)
        
        -- reactions from model, allow only ECs from model
        WHERE Reaction_ID_MetaCyc IN (%s)
        AND e.ec_number IN (%s)
        """ % (reactionString, ecString))
        
        result = cursor.fetchall()
        
        if len(result)<=0:
            print("No reaction could be mapped.")
            return
        
        self.drawReactionNodes(list(map(str,flatten(result))))
        
        # find untranslated reaction ids
        cursor.execute("""
        SELECT r2db.Reaction_ID_MetaCyc
        FROM reaction2kegg_metacyc_brenda AS r2db
        """)
        
        result = cursor.fetchall()
        allDBids = set(flatten(result))
        missingIDs = set(rid for rid, _ in reactionEClist).difference(allDBids)
#        pprint(missingIDs)
        
        reactionIDtoECdict = dict((r,e) for r,e in reactionEClist if e is not None and r in missingIDs)
        
#        pprint(reactionIDtoECdict)
        self.drawEnzymes(list(reactionIDtoECdict.values()))
        
        # TODO return list of drawn reactions

    def drawRemainingReactions(self, alpha="30"): # alpha unused
        
        connection = self.__getConnection()
        cursor = connection.cursor()

        if len(self.graph)<=0:
            print("The graph is empty.")
            return

        nodeString = listToString(self.graph.nodes())
        
        # add all edges and nodes to the graph
        
        cursor.execute("""
        SELECT DISTINCT node_id
        FROM enzyme2node
        WHERE node_id NOT IN (%s)
        """ % nodeString)
        
        result = cursor.fetchall()
        
        # TODO set as parameter(object?)
        self.reactionColor = "#ffcd0030"
        self.substanceColor = "#89a40030"
        self.edgeColor = "#00000030"
        
        self.drawReactionNodes(flatten(result))

    def drawPathwayLabels(self):

        connection = self.__getConnection()
        cursor = connection.cursor()

        # draw pathway frame nodes

#        cursor.execute("""
#            SELECT MIN(node_id) AS node_id, AVG(x_pos) AS x, AVG(y_pos) AS y,
#            (MAX(x_pos)-MIN(x_pos)) AS width, (MAX(y_pos)-MIN(y_pos)) AS height
#            FROM `node`
#            WHERE `node_type` LIKE 'pathwayFrame'
#            GROUP BY `dominant_pathway_id`
#            """)
#        result = cursor.fetchall()
#        
#        for nodeID, x, y, width, height in result:
#            self.graph.add_node(nodeID)
#            self.graph.node[nodeID]["label"]=""
#            self.graph.node[nodeID]["graphics"]= \
#                {"x":x,
#                 "y":y,
#                 "w":width,
#                 "h":height,
#                 "hasFill":"0",
#                 "outlineStyle": "\"dashed\"",
#                 "type":"\"rectangle\"",
#                 }

        # draw pathway label nodes
        # TODO increase font size
        # TODO color

        cursor.execute("""
            SELECT node_id, node_label, x_pos, y_pos, width, height
            FROM node
            WHERE `node_type` LIKE 'Pathway_title'
            """)
            
        result = cursor.fetchall()
        
        for nodeID, label, x, y, width, height in result:
            self.graph.add_node(nodeID)
            self.graph.node[nodeID]["label"]=label
            self.graph.node[nodeID]["graphics"]= \
                {"x":x,
                 "y":float(y)-0.5*30., # move labels
                 "w":max(30.0, len(label)*9.), # TODO magic numbers
                 "h":30.0,
                 "type":self.substanceShape,
#                 "fill":"\""+self.reactionColor+"\""
                 }

    def setNodeGraphicsAttribute(self, attribute, valueList, attributeDict):
        """
        Set all graphics attributes in attributeDict of all nodes in the graph
        whose value of the given attribute is in the valueList
        
        @param attribute: The nodes will be filtered by this attribute
        @type attribute: string
        @param valueList: List of values to match with the attribute
        @type valueList: list
        @param attributeDict: Dict of GML graphics attributes
        @type attributeDict: dict
        
        returns a subset of valueList with attributes which could _not_ be
        found in any node of the graph.
        """
        # TODO example

        foundAttributes = set()
        
        for node in self.graph.nodes_iter():

            if self.graph.node[node][attribute] in valueList:

                if self.verbosityLevel>1:
                    print(node, "with attribute", \
                          self.graph.node[node][attribute], "found")
                
                foundAttributes.add(self.graph.node[node][attribute])              
            
                for key, value in attributeDict.items():
                    self.graph.node[node]["graphics"][key]=value

        if self.verbosityLevel>0:
            print(len(foundAttributes), "of", len(set(valueList)), "distinct", \
                  "attributes found")

        return set(valueList).difference(foundAttributes)


def main(argv=None):

    if argv is None:
        argv = sys.argv

    server = "fileserver.bioinfo.nat.tu-bs.de"
    user = "rre"
    password = open("pw.txt").readline().strip()
    database = "metabolic_pathways"
    verbosityLevel = 1

#    mymap = MetabolicMap(server, user, password, database, verbosityLevel)
#    mymap.drawEnzymes(["6.4.1.1", "1.1.1.35"])
#    mymap.writeToGML("/home/rre/Desktop/test.gml")
#    mymap.writeToDot("/home/rre/Desktop/test.dot")
#    mymap.writeToGraphML("/home/rre/Desktop/test.graphml")
#    mymap.closeConnection()
    
    connection = MySQLdb.connect(server, user,password,"rre_roseobacter")
    cursor = connection.cursor()
    
    megateriumMap = MetabolicMap(server, user, password, database, verbosityLevel)
    
#    cursor.execute("""
#    SELECT DISTINCT Bmegaterium__model.ec
#    FROM Bmegaterium__model
#    INNER JOIN Bmegaterium__scenario
#    ON (Bmegaterium__model.reaction_id=Bmegaterium__scenario.reaction_id)""")
#    result = cursor.fetchall()
#    ecList = flatten(result)
#    megateriumMap.drawEnzymes(ecList)

#    cursor.execute("""
#    SELECT DISTINCT reaction_id
#    FROM Bmegaterium__scenario
#    """)
#    result = cursor.fetchall()
#    reactionIDlist = flatten(result)
#    megateriumMap.drawReactions(reactionIDlist)

    cursor.execute("""
    SELECT reaction_id, ec
    FROM Bmegaterium__table_BiomassGFP
    """)
    reactionEClist = cursor.fetchall()
    megateriumMap.drawReactionEnzyme(reactionEClist)
    megateriumMap.drawRemainingReactions()
    megateriumMap.drawPathwayLabels()

    # TODO metabolite liste einlesen
#    metabList=list(open("/home/rre/Desktop/Megaterium/BmegateriumT2.csv"))
#    metabList=[m.strip() for m in metabList]
#    notFound = megateriumMap.setNodeGraphicsAttribute("label",metabList,{"fill":"\"#CC0000\""})
#    print notFound
    

    megateriumMap.writeToGML("/home/rre/Desktop/Bmega.gml")
    
    megateriumMap.closeConnection()
    connection.close()
    
#    dinoMap = MetabolicMap(server, user, password, database, verbosityLevel)    
#    cursor.execute("""
#    SELECT reaction_id, ec
#    FROM `Dshibae__table_beta-D-glucose`
#    """)
#    reactionEClist = cursor.fetchall()
#    dinoMap.drawReactionEnzyme(reactionEClist)
#    dinoMap.drawPathwayLabels()
#    dinoMap.drawRemainingReactions()
#    dinoMap.writeToGML("/home/rre/Desktop/dino-glucose.gml")
#    dinoMap.closeConnection()
#    
#    dinoMap = MetabolicMap(server, user, password, database, verbosityLevel)
#    cursor.execute("""
#    SELECT reaction_id, ec
#    FROM Dshibae__model
#    """)
#    reactionEClist = cursor.fetchall()
#    dinoMap.drawReactionEnzyme(reactionEClist)
#    dinoMap.drawPathwayLabels()
#    dinoMap.drawRemainingReactions()
#    dinoMap.writeToGML("/home/rre/Desktop/dino-all.gml")
#    dinoMap.closeConnection()
    
#    phaeoMap = MetabolicMap(server, user, password, database, verbosityLevel)
#    cursor.execute("""
#    SELECT reaction_id, ec
#--    FROM `Pgallaeciensis__succinate_L-valine`
#    FROM `Pgallaeciensis__succinate_beta-D-glucose`
#    """)
#    reactionEClist = cursor.fetchall()
#    phaeoMap.drawReactionEnzyme(reactionEClist)
#    phaeoMap.drawPathwayLabels()
#    phaeoMap.drawRemainingReactions()
#    phaeoMap.writeToGML("/home/rre/Desktop/phaeo-glucose.gml")
#    phaeoMap.closeConnection()
    
    return 0

if __name__ == "__main__":

    sys.exit(main())
