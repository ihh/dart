import sys
class SExpr(list):
    # A fairly wimpy SExpr class
    def make_all_SExprs(self):
        for i,item in enumerate(self):
            if type(item) == list:
                toInsert = self.pop(i)
                self.insert(i,SExpr(toInsert))
                self[i].make_all_SExprs()

    def has_tag(self):
        for item in self:
            if type(item) != list:
                return True
        return False
    def get_tag(self):
        for item in self:
            if type(item) != list:
                return item
    def is_tag_value_pair(self):
        return len(self) == 2
    def find_all_with_tag(self, query):
        out = []
        for item in self:
            if type(item)==str:
                continue
            if item.has_tag():
                if item.get_tag() == query:
                    out.append(item)
        if not out:
            sys.stderr.write("Warning: no children matched the tag " +query+'\n')
        return out
    def get_value_for_tag(self, query):
        for item in self:
            if type(item)==str:
                continue
            if item.has_tag():
                if item.get_tag() == query:
                    return item[1]
        sys.stderr.write("Warning: no entry found for tag " + query +'\n')

            
    

        
        
