package beast.app.beastlier.beauti;

import beast.app.beauti.Beauti;
import beast.app.util.Utils;
import beast.core.util.Log;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;


public class GuessPatternDialogBEASTLIER extends JPanel {
    private static final long serialVersionUID = 1L;

    public final String helpText;

    public enum Status {
        canceled, pattern, trait
    };

    public String trait = null;

    public String getTrait() {
        return trait;
    }

    Component m_parent;
    JPanel guessPanel;
    ButtonGroup group;
    JRadioButton useEverything = new JRadioButton("use everything");
    JRadioButton isSplitOnChar = new JRadioButton("split on character");
    JRadioButton useRegexp = new JRadioButton("use regular expression");
    JRadioButton readFromFile = new JRadioButton("read from");

    int m_location = 0;
    int m_splitlocation = 0;
    int m_textFileDelimiter = 0;
    String m_sDelimiter = ".";
    JTextField textRegExp;
    JComboBox<String> combo;
    JComboBox<String> combo_1;
    JComboBox<String> combo_2;
    String pattern;
    String textFileDelimiter;

    public String getPattern() {
        return pattern;
    }

    private JTextField txtFile;
    private JTextField textSplitChar;
    private JTextField textSplitChar2;
    private JTextField textAddValue;
    private JTextField textUnlessLessThan;
    private JTextField textThenAdd;
    JCheckBox chckbxAddFixedValue;
    JCheckBox chckbxUnlessLessThan;
    JLabel lblThenAdd;
    JLabel lblAndTakeGroups;
    JLabel lblFile;
    JButton btnBrowse;
    private JSeparator separator_2;
    private JSeparator separator_3;
    private JSeparator separator_4;
    private JSeparator separator_5;

    private String indexColumnHeader;
    private String valueColumnHeader;

    public GuessPatternDialogBEASTLIER(Component parent, String pattern, String indexColumnHeader,
                                       String valueColumnHeader) {
        helpText = "<html>The loaded file can be in tab-separated or comma-separated format. The first row <br>" +
                "must contain the column headers, which must include <font color=red>"+indexColumnHeader+"</font>"
                + "and <font color=red>"+valueColumnHeader+"</font> <br>(in any order). Other columns will be ignored.";
        this.indexColumnHeader = indexColumnHeader;
        this.valueColumnHeader = valueColumnHeader;

        m_parent = parent;
        this.pattern = pattern;
        guessPanel = new JPanel();
        GridBagLayout gbl_guessPanel = new GridBagLayout();
        gbl_guessPanel.rowHeights = new int[]{0, 0, 0, 20, 0, 0, 20, 0, 0, 20, 0, 29, 0, 0, 0, 0};
        gbl_guessPanel.columnWidths = new int[]{0, 0, 0, 0, 0, 0};
        gbl_guessPanel.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        gbl_guessPanel.columnWeights = new double[] { 1.0, 1.0, 1.0, 0.0, 1.0, 0.0 };
        guessPanel.setLayout(gbl_guessPanel);

        group = new ButtonGroup();
        group.add(useEverything);
        group.add(isSplitOnChar);
        group.add(useRegexp);
        group.add(readFromFile);
        group.setSelected(useEverything.getModel(), true);
        useEverything.addActionListener(e -> {
            updateFields();
        });
        useEverything.setName(useEverything.getText());
        isSplitOnChar.addActionListener(e -> {
            updateFields();
        });
        isSplitOnChar.setName(isSplitOnChar.getText());
        useRegexp.addActionListener(e -> {
            updateFields();
        });
        useRegexp.setName(useRegexp.getText());
        readFromFile.addActionListener(e -> {
            updateFields();
        });
        readFromFile.setName(readFromFile.getText());

        createDelimiterBox(useEverything);
        createSplitBox(isSplitOnChar);
        createRegExpBox(useRegexp);

        textRegExp = new JTextField();
        textRegExp.setText(pattern);
        textRegExp.setColumns(10);
        textRegExp.setToolTipText("Enter regular expression to match taxa");
        textRegExp.setMaximumSize(new Dimension(1024, 25));
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.insets = new Insets(0, 0, 5, 5);
        gbc2.anchor = GridBagConstraints.WEST;
        gbc2.gridwidth = 4;
        gbc2.gridx = 1;
        gbc2.gridy = 7;
        guessPanel.add(textRegExp, gbc2);
        textRegExp.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void removeUpdate(DocumentEvent e) {
                useRegexp.setSelected(true);
            }

            @Override
            public void insertUpdate(DocumentEvent e) {
                useRegexp.setSelected(true);
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                useRegexp.setSelected(true);
            }
        });

        separator_4 = new JSeparator();
        separator_4.setPreferredSize(new Dimension(5,1));
        GridBagConstraints gbc_separator_4 = new GridBagConstraints();
        gbc_separator_4.gridwidth = 5;
        gbc_separator_4.insets = new Insets(5, 0, 15, 5);
        gbc_separator_4.gridx = 0;
        gbc_separator_4.gridy = 8;
        gbc_separator_4.fill = GridBagConstraints.HORIZONTAL;
        guessPanel.add(separator_4, gbc_separator_4);

        GridBagConstraints gbc_rdbtnReadFromFile = new GridBagConstraints();
        gbc_rdbtnReadFromFile.anchor = GridBagConstraints.WEST;
        gbc_rdbtnReadFromFile.insets = new Insets(0, 0, 5, 5);
        gbc_rdbtnReadFromFile.gridx = 0;
        gbc_rdbtnReadFromFile.gridy = 10;
        guessPanel.add(readFromFile, gbc_rdbtnReadFromFile);

        combo_2 = new JComboBox<>(new String[]{".txt", ".csv"});
        combo_2.setName("fileType");
        GridBagConstraints gbc_format = new GridBagConstraints();
        gbc_format.gridwidth = 1;
        gbc_format.insets = new Insets(0, 0, 5, 5);
        gbc_format.fill = GridBagConstraints.HORIZONTAL;
        gbc_format.gridx = 1;
        gbc_format.gridy = 10;

        combo_2.addActionListener(e -> {
            @SuppressWarnings("unchecked")
            JComboBox<String> combo = (JComboBox<String>) e.getSource();
            m_textFileDelimiter = combo.getSelectedIndex();
            readFromFile.setSelected(true);
            updateFields();
        });
        guessPanel.add(combo_2, gbc_format);

        lblFile = new JLabel("file");
        GridBagConstraints gbc_lblFile = new GridBagConstraints();
        gbc_lblFile.gridwidth = 2;
        gbc_lblFile.insets = new Insets(0, 0, 5, 5);
        gbc_lblFile.gridx = 2;
        gbc_lblFile.gridy = 10;
        guessPanel.add(lblFile, gbc_lblFile);

        btnBrowse = new JButton("Browse");
        btnBrowse.addActionListener(e -> {
            File file = Utils.getLoadFile("Load trait from file", new File(Beauti.g_sDir), "Select trait file", "dat",
                    "txt","csv");
            if (file != null) {
                txtFile.setText(file.getPath());
                readFromFile.setSelected(true);
                updateFields();
            }
        });

        txtFile = new JTextField();
        txtFile.setText("file");
        GridBagConstraints gbc_txtFile = new GridBagConstraints();
        gbc_txtFile.gridwidth = 2;
        gbc_txtFile.insets = new Insets(0, 0, 5, 5);
        gbc_txtFile.fill = GridBagConstraints.HORIZONTAL;
        gbc_txtFile.gridx = 3;
        gbc_txtFile.gridy = 10;
        guessPanel.add(txtFile, gbc_txtFile);
        txtFile.setColumns(10);
        GridBagConstraints gbc_btnReadFromFile = new GridBagConstraints();
        gbc_btnReadFromFile.insets = new Insets(0, 0, 5, 5);
        gbc_btnReadFromFile.gridx = 5;
        gbc_btnReadFromFile.gridy = 10;
        guessPanel.add(btnBrowse, gbc_btnReadFromFile);

        JButton btnHelp = new JButton("?");
        btnHelp.setToolTipText("Show format of trait file");
        btnHelp.addActionListener(e -> {
            JOptionPane.showMessageDialog(m_parent, helpText);
        });
        GridBagConstraints gbc_btnHelp = new GridBagConstraints();
        gbc_btnHelp.insets = new Insets(0, 0, 5, 5);
        gbc_btnHelp.gridx = 6;
        gbc_btnHelp.gridy = 10;
        guessPanel.add(btnHelp, gbc_btnHelp);


        chckbxAddFixedValue = new JCheckBox("Add fixed value");
        chckbxAddFixedValue.setName("Add fixed value");
        chckbxAddFixedValue.addActionListener(e -> {
            updateFields();
        });

        separator_5 = new JSeparator();
        separator_5.setPreferredSize(new Dimension(5,1));
        GridBagConstraints gbc_separator_5 = new GridBagConstraints();
        gbc_separator_5.gridwidth = 5;
        gbc_separator_5.insets = new Insets(5, 0, 15, 5);
        gbc_separator_5.gridx = 0;
        gbc_separator_5.gridy = 12;
        gbc_separator_5.fill = GridBagConstraints.HORIZONTAL;
        guessPanel.add(separator_5, gbc_separator_5);
        GridBagConstraints gbc_chckbxAddFixedValue = new GridBagConstraints();
        gbc_chckbxAddFixedValue.anchor = GridBagConstraints.WEST;
        gbc_chckbxAddFixedValue.insets = new Insets(0, 0, 5, 5);
        gbc_chckbxAddFixedValue.gridx = 0;
        gbc_chckbxAddFixedValue.gridy = 13;
        guessPanel.add(chckbxAddFixedValue, gbc_chckbxAddFixedValue);

        textAddValue = new JTextField("1900");
        GridBagConstraints gbc_textField = new GridBagConstraints();
        gbc_textField.gridwidth = 2;
        gbc_textField.insets = new Insets(0, 0, 5, 5);
        gbc_textField.fill = GridBagConstraints.HORIZONTAL;
        gbc_textField.gridx = 1;
        gbc_textField.gridy = 13;
        guessPanel.add(textAddValue, gbc_textField);
        textAddValue.setColumns(10);

        chckbxUnlessLessThan = new JCheckBox("Unless less than...");
        chckbxUnlessLessThan.setName("Unless less than");
        chckbxUnlessLessThan.addActionListener(e -> {
            updateFields();
        });
        GridBagConstraints gbc_chckbxUnlessLargerThan = new GridBagConstraints();
        gbc_chckbxUnlessLargerThan.anchor = GridBagConstraints.WEST;
        gbc_chckbxUnlessLargerThan.insets = new Insets(0, 0, 5, 5);
        gbc_chckbxUnlessLargerThan.gridx = 0;
        gbc_chckbxUnlessLargerThan.gridy = 14;
        guessPanel.add(chckbxUnlessLessThan, gbc_chckbxUnlessLargerThan);

        textUnlessLessThan = new JTextField("13");
        GridBagConstraints gbc_textField_1 = new GridBagConstraints();
        gbc_textField_1.gridwidth = 2;
        gbc_textField_1.insets = new Insets(0, 0, 5, 5);
        gbc_textField_1.fill = GridBagConstraints.HORIZONTAL;
        gbc_textField_1.gridx = 1;
        gbc_textField_1.gridy = 14;
        guessPanel.add(textUnlessLessThan, gbc_textField_1);
        textUnlessLessThan.setColumns(10);

        lblThenAdd = new JLabel("...then add");
        GridBagConstraints gbc_lblThenAdd = new GridBagConstraints();
        gbc_lblThenAdd.anchor = GridBagConstraints.EAST;
        gbc_lblThenAdd.insets = new Insets(0, 0, 0, 5);
        gbc_lblThenAdd.gridx = 0;
        gbc_lblThenAdd.gridy = 15;
        guessPanel.add(lblThenAdd, gbc_lblThenAdd);

        textThenAdd = new JTextField("2000");
        GridBagConstraints gbc_textField_2 = new GridBagConstraints();
        gbc_textField_2.gridwidth = 2;
        gbc_textField_2.insets = new Insets(0, 0, 0, 5);
        gbc_textField_2.fill = GridBagConstraints.HORIZONTAL;
        gbc_textField_2.gridx = 1;
        gbc_textField_2.gridy = 15;
        guessPanel.add(textThenAdd, gbc_textField_2);
        textThenAdd.setColumns(10);


        chckbxAddFixedValue.setVisible(false);
        textAddValue.setVisible(false);
        chckbxUnlessLessThan.setVisible(false);
        lblThenAdd.setVisible(false);
        chckbxUnlessLessThan.setVisible(false);
        textUnlessLessThan.setVisible(false);
        textThenAdd.setVisible(false);
    }

    public void allowAddingValues() {
        chckbxAddFixedValue.setVisible(true);
        textAddValue.setVisible(true);
        chckbxUnlessLessThan.setVisible(true);
        lblThenAdd.setVisible(true);
        chckbxUnlessLessThan.setVisible(true);
        textUnlessLessThan.setVisible(true);
        textThenAdd.setVisible(true);
    }

    protected void updateFields() {
        if (chckbxAddFixedValue.isSelected()) {
            textAddValue.setEnabled(true);
            chckbxUnlessLessThan.setEnabled(true);
            lblThenAdd.setEnabled(true);
            if (chckbxUnlessLessThan.isSelected()) {
                textUnlessLessThan.setEnabled(true);
                textThenAdd.setEnabled(true);
            } else {
                textUnlessLessThan.setEnabled(false);
                textThenAdd.setEnabled(false);
            }
        } else {
            textAddValue.setEnabled(false);
            chckbxUnlessLessThan.setEnabled(false);
            lblThenAdd.setEnabled(false);
            textUnlessLessThan.setEnabled(false);
            textThenAdd.setEnabled(false);
        }

        txtFile.setEnabled(false);
        textSplitChar.setEnabled(false);
        textSplitChar2.setEnabled(false);
        textRegExp.setEnabled(false);
        combo.setEnabled(false);
        combo_1.setEnabled(false);
        lblAndTakeGroups.setEnabled(false);
        btnBrowse.setEnabled(false);
        if (useEverything.isSelected()) {
            textSplitChar.setEnabled(true);
            combo.setEnabled(true);
        }
        if (isSplitOnChar.isSelected()) {
            textSplitChar2.setEnabled(true);
            combo_1.setEnabled(true);
            lblAndTakeGroups.setEnabled(true);
        }
        if (useRegexp.isSelected()) {
            textRegExp.setEnabled(true);
        }
        if (readFromFile.isSelected()) {
            btnBrowse.setEnabled(true);
            txtFile.setEnabled(true);
        }
    }

    private void createDelimiterBox(JRadioButton b) {
        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(0, 0, 5, 5);
        gbc.anchor = GridBagConstraints.WEST;
        gbc.gridx = 0;
        gbc.gridy = 1;
        guessPanel.add(b, gbc);

        combo = new JComboBox<>(new String[] { "after first", "after last", "before first", "before last" });
        combo.setName("delimiterCombo");
        GridBagConstraints gbc2 = new GridBagConstraints();
        gbc2.anchor = GridBagConstraints.WEST;
        gbc2.gridwidth = 2;
        gbc2.insets = new Insets(0, 0, 5, 5);
        gbc2.gridx = 1;
        gbc2.gridy = 1;
        guessPanel.add(combo, gbc2);
        combo.addActionListener(e -> {
            @SuppressWarnings("unchecked")
            JComboBox<String> combo = (JComboBox<String>) e.getSource();
            m_location = combo.getSelectedIndex();
            useEverything.setSelected(true);
            updateFields();
        });
    }

    private void createSplitBox(JRadioButton b) {

        textSplitChar = new JTextField("_");
        textSplitChar.setName("SplitChar");
        GridBagConstraints gbc_textField = new GridBagConstraints();
        gbc_textField.anchor = GridBagConstraints.WEST;
        gbc_textField.insets = new Insets(0, 0, 5, 5);
        gbc_textField.gridx = 3;
        gbc_textField.gridy = 1;
        guessPanel.add(textSplitChar, gbc_textField);
        textSplitChar.setColumns(2);

        separator_2 = new JSeparator();
        separator_2.setPreferredSize(new Dimension(5,1));
        GridBagConstraints gbc_separator_2 = new GridBagConstraints();
        gbc_separator_2.gridwidth = 5;
        gbc_separator_2.insets = new Insets(5, 0, 15, 5);
        gbc_separator_2.gridx = 0;
        gbc_separator_2.gridy = 2;
        gbc_separator_2.fill = GridBagConstraints.HORIZONTAL;
        guessPanel.add(separator_2, gbc_separator_2);

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(0, 0, 5, 5);
        gbc.anchor = GridBagConstraints.WEST;
        gbc.gridx = 0;
        gbc.gridy = 4;
        guessPanel.add(b, gbc);
    }

    public void createRegExpBox(JRadioButton b) {

        textSplitChar2 = new JTextField("_");
        textSplitChar2.setName("SplitChar2");
        GridBagConstraints gbc_textField_1 = new GridBagConstraints();
        gbc_textField_1.anchor = GridBagConstraints.WEST;
        gbc_textField_1.insets = new Insets(0, 0, 5, 5);
        gbc_textField_1.gridx = 1;
        gbc_textField_1.gridy = 4;
        guessPanel.add(textSplitChar2, gbc_textField_1);
        textSplitChar2.setColumns(2);

        lblAndTakeGroups = new JLabel("and take group(s):");
        GridBagConstraints gbc_lblAndTakeGroups = new GridBagConstraints();
        gbc_lblAndTakeGroups.gridwidth = 2;
        gbc_lblAndTakeGroups.insets = new Insets(0, 0, 5, 5);
        gbc_lblAndTakeGroups.gridx = 2;
        gbc_lblAndTakeGroups.gridy = 4;
        guessPanel.add(lblAndTakeGroups, gbc_lblAndTakeGroups);

        combo_1 = new JComboBox<>(new String[] { "1", "2", "3", "4", "1-2", "2-3", "3-4", "1-3", "2-4" });
        combo_1.setName("splitCombo");
        GridBagConstraints gbc_combo_1 = new GridBagConstraints();
        gbc_combo_1.anchor = GridBagConstraints.WEST;
        gbc_combo_1.insets = new Insets(0, 0, 5, 5);
        gbc_combo_1.gridx = 4;
        gbc_combo_1.gridy = 4;
        guessPanel.add(combo_1, gbc_combo_1);
        combo_1.addActionListener(e -> {
            @SuppressWarnings("unchecked")
            JComboBox<String> combo = (JComboBox<String>) e.getSource();
            m_splitlocation = combo.getSelectedIndex();
            isSplitOnChar.setSelected(true);
            updateFields();
        });

        separator_3 = new JSeparator();
        separator_3.setPreferredSize(new Dimension(5,1));
        GridBagConstraints gbc_separator_3 = new GridBagConstraints();
        gbc_separator_3.gridwidth = 5;
        gbc_separator_3.insets = new Insets(5, 0, 15, 5);
        gbc_separator_3.gridx = 0;
        gbc_separator_3.gridy = 5;
        gbc_separator_3.fill = GridBagConstraints.HORIZONTAL;
        guessPanel.add(separator_3, gbc_separator_3);

        GridBagConstraints gbc = new GridBagConstraints();
        gbc.insets = new Insets(0, 0, 5, 5);
        gbc.anchor = GridBagConstraints.WEST;
        gbc.gridx = 0;
        gbc.gridy = 7;
        guessPanel.add(b, gbc);
    }

    public Status showDialog(String title) {

        JOptionPane optionPane = new JOptionPane(guessPanel, JOptionPane.PLAIN_MESSAGE, JOptionPane.OK_CANCEL_OPTION,
                null, new String[] { "Cancel", "OK" }, "OK");
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(m_parent, title);
        dialog.setName("GuessTaxonSets");
        // dialog.setResizable(true);
        dialog.pack();
        updateFields();
        dialog.setVisible(true);

        if (optionPane.getValue() == null || !optionPane.getValue().equals("OK")) {
            return Status.canceled;
        }

        if (useEverything.getModel() == group.getSelection()) {
            String delimiter = normalise(textSplitChar.getText());
            switch (m_location) {
                case 0: // "after first",
                    pattern = "^[^" + delimiter + "]+" + delimiter + "(.*)$";
                    break;
                case 1: // "after last",
                    pattern = "^.*" + delimiter + "(.*)$";
                    break;
                case 2: // "before first",
                    pattern = "^([^" + delimiter + "]+)" + delimiter + ".*$";
                    break;
                case 3: // "before last"
                    pattern = "^(.*)" + delimiter + ".*$";
                    break;
            }
        }
        if (isSplitOnChar.getModel() == group.getSelection()) {
            String delimiter = normalise(textSplitChar2.getText());
            switch (m_splitlocation) {
                case 0: // "1"
                    pattern = "^([^" + delimiter + "]+)" + ".*$";
                    break;
                case 1: // "2"
                    pattern = "^[^" + delimiter + "]+" + delimiter + "([^" + delimiter + "]+)" + ".*$";
                    break;
                case 2: // "3"
                    pattern = "^[^" + delimiter + "]+" + delimiter + "[^" + delimiter + "]+" + delimiter + "([^"
                            + delimiter + "]+)" + ".*$";
                    break;
                case 3: // "4"
                    pattern = "^[^" + delimiter + "]+" + delimiter + "[^" + delimiter + "]+" + delimiter + "[^"
                            + delimiter + "]+" + delimiter + "([^" + delimiter + "]+)" + ".*$";
                    break;
                case 4: // "1-2"
                    pattern = "^([^" + delimiter + "]+" + delimiter + "[^" + delimiter + "]+)" + ".*$";
                    break;
                case 5: // "2-3"
                    pattern = "^[^" + delimiter + "]+" + delimiter + "([^" + delimiter + "]+" + delimiter + "[^"
                            + delimiter + "]+)" + ".*$";
                    break;
                case 6: // "3-4"
                    pattern = "^[^" + delimiter + "]+" + delimiter + "[^" + delimiter + "]+" + delimiter + "([^"
                            + delimiter + "]+" + delimiter + "[^" + delimiter + "]+)" + ".*$";
                    break;
                case 7: // "1-3"
                    pattern = "^([^" + delimiter + "]+" + delimiter + "[^" + delimiter + "]+" + delimiter + "[^"
                            + delimiter + "]+)" + ".*$";
                    break;
                case 8: // "2-4"
                    pattern = "^[^" + delimiter + "]+" + delimiter + "([^" + delimiter + "]+" + delimiter + "[^"
                            + delimiter + "]+" + delimiter + "[^" + delimiter + "]+)" + ".*$";
            }
        }
        if (useRegexp.getModel() == group.getSelection()) {
            pattern = textRegExp.getText();
        }
        if (readFromFile.getModel() == group.getSelection()) {
            switch (m_textFileDelimiter){
                case 0:
                    textFileDelimiter = "\t";
                    break;
                case 1:
                    textFileDelimiter = ",";
            }
            try {
                BufferedReader fin = new BufferedReader(new FileReader(txtFile.getText()));
                StringBuffer buf = new StringBuffer();
                //the header
                String header = fin.readLine();
                String[] splitHeader = header.split(textFileDelimiter);
                int indexCol = -1;
                int valueCol = -1;

                for(int i=0; i<splitHeader.length; i++){
                    if(splitHeader[i].equals(indexColumnHeader)){
                        indexCol = i;
                    }
                    if(splitHeader[i].equals(valueColumnHeader)){
                        valueCol = i;
                    }
                }

                while (fin.ready()) {
                    String line = fin.readLine();
                    String[] splitLine = line.split(textFileDelimiter);
                    String str = splitLine[indexCol]+"="+splitLine[valueCol]+",";
                    // only add entries that are non-empty
                    if (str.indexOf("=") > 0 && !str.matches("^\\s+=.*$")) {
                        buf.append(str);
                    }
                }
                fin.close();
                trait = buf.toString().trim();
                while (trait.endsWith(",")) {
                    trait = trait.substring(0, trait.length() - 1).trim();
                }

                if (indexCol == -1) {
                    JOptionPane.showMessageDialog(m_parent, "Could not find the required header "+indexColumnHeader+
                                    " in this file");
                }
                if (valueCol == -1) {
                    JOptionPane.showMessageDialog(m_parent, "Could not find the required header "+valueColumnHeader+
                            " in this file");
                }
            } catch (Exception e) {
                JOptionPane.showMessageDialog(m_parent, "Loading trait from file failed:" + e.getMessage());
                return Status.canceled;
            }
            return Status.trait;
        }

        // sanity check
        try {
            pattern.matches(pattern);
        } catch (PatternSyntaxException e) {
            JOptionPane.showMessageDialog(this, "This is not a valid regular expression");
            return Status.canceled;
        }

        if (optionPane.getValue() != null && optionPane.getValue().equals("OK")) {
            Log.warning.println("Pattern = " + pattern);
            return Status.pattern;
        } else {
            return Status.canceled;
        }
    }

    /**
     * Converts the first character of delimiter into a substring suitable for
     * inclusion in a regexp. This is done by expressing the character as an
     * octal escape.
     *
     * @param delimiter first character of this string to be used as delimiter
     * @return escaped octal representation of character
     */
    private String normalise(String delimiter) {

        if (delimiter.length() == 0) {
            return ".";
        }

        return String.format("\\0%o", (int)delimiter.charAt(0));
    }

    public String match(String s) {
        Pattern _pattern = Pattern.compile(pattern);
        Matcher matcher = _pattern.matcher(s);
        if (matcher.find()) {
            String match = matcher.group(1);
            if (chckbxAddFixedValue.isSelected()) {
                try {
                    Double value = Double.parseDouble(match);
                    Double addValue = Double.parseDouble(textAddValue.getText());
                    if (chckbxUnlessLessThan.isSelected()) {
                        Double threshold = Double.parseDouble(textUnlessLessThan.getText());
                        Double addValue2 = Double.parseDouble(textThenAdd.getText());
                        if (value < threshold) {
                            value += addValue2;
                        } else {
                            value += addValue;
                        }
                    } else {
                        value += addValue;
                    }
                    return value + "";
                } catch (Exception e) {
                    // ignore
                }
            }
            return match;
        }
        return null;
    }
}
