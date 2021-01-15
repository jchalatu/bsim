package PhysModBsim;

import java.util.ArrayList;
import java.util.List;

public class Node<T> {
    private T data;
    private List<Node<T>> children;
    private Node<T> parent;


    public Node(T data) {
        this.data = data;
        children = new ArrayList<Node<T>>();
        parent = null;
    }

    public Node(T data, Node<T> parent) {
        this.data = data;
        children = new ArrayList<Node<T>>();
        this.parent = parent;
        parent.children.add(this);
    }

    /*
    public void addChild(Node<T> child) {
        children.add(child);
        child.parent = this;
    }

    public void addChild(T data) {
        Node<T> child = new Node<T>(data);
        children.add(child);
        child.parent = this;
    }
    */

    public T getData() {
        return data;
    }

    public void setData(T data) {
        this.data = data;
    }

    public List<Node<T>> getChildren() {
        return children;
    }

    public boolean isRoot() {
        return (this.parent == null);
    }

    public boolean isLeaf() {
        return this.children.size()  == 0;
    }
}
